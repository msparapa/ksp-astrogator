import time
import beluga
import krpc
import logging
from math import pi, sqrt
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

test_run = False
countdown = 3

conn = krpc.connect(name='ksp-astrogator')

if not test_run:
    vessel = conn.space_center.active_vessel
    print('Connected with ' + vessel.name)


ocp = beluga.OCP('goddard')

# Define independent variables
ocp.independent('t', 's')

# Define quantities used in the problem
ocp.quantity('drag', '1/2 * d * v**2 * rho0 * exp(-h/H) * pi * (A/2)**2')

# Define equations of motion
ocp.state('h', 'v', 'm') \
   .state('v', '(thrust*(thrust_max_sea*exp(-h/H) + thrust_max_vac*(1-exp(-h/H))) - drag)/m - g', 'm/s') \
   .state('m', '-thrust*c', 'kg')

# Define controls
ocp.control('thrust', '1')

# Define constants' numerical values
g_0 = 9.81
h_0 = 0
v_0 = 0
m_0 = 6.275e3
rho0 = 0.00
A = 1.4

c = 68.51
H = 5600
tar_m_f = 2.275e3
d = 0.2
drag_max = 17500
thrust_max_sea = 167.97e3
thrust_max_vac = 215.00e3

meter = H
second = 5
kilogram = m_0

# Define constants
ocp.constant('g', g_0/meter*second**2, '1')  # Gravity at surface
ocp.constant('H', H/meter, '1')      # Constant for height

ocp.constant('c', c/kilogram*second, '1')          # Thrust to fuel ratio
ocp.constant('d', d, '1')      # Drag scaling
ocp.constant('thrust_max_sea', thrust_max_sea/kilogram/meter*second**2, '1')
ocp.constant('thrust_max_vac', thrust_max_vac/kilogram/meter*second**2, '1')
ocp.constant('rho0', rho0/kilogram*meter**3, '1')
ocp.constant('drag_max', drag_max/kilogram/meter*second**2, '1')

# Define constants for BCs
ocp.constant('h_0', h_0/meter, '1')
ocp.constant('v_0', v_0/meter*second, '1')
ocp.constant('m_0', m_0/m_0, '1')
ocp.constant('A', A/meter, '1')

ocp.constant('v_f', v_0/H*second, '1')
ocp.constant('m_f', tar_m_f/kilogram, '1')

# Define smoothing constant
ocp.constant('eps', 0.1, '1')
ocp.constant('eps2', 0.005, '1')

# Define costs
ocp.terminal_cost('-h', '1')
# ocp.path_cost('-eps*cos(u)', '1')

# Define constraints
ocp.constraints() \
    .initial('h - h_0', '1') \
    .initial('v - v_0', '1') \
    .initial('m - m_0', '1') \
    .path('thrust', '1', lower=0, upper=1, activator='eps', method='epstrig') \
    .path('1/2 * d * v**2 * rho0 * exp(-h/H) * pi * (A/2)**2', '1', lower=-1, upper='drag_max', activator='eps2',
          method='utm') \
    .terminal('v - v_f', '1') \
    .terminal('m - m_f', '1')


ocp.scale(m=1, s=1, kg=1, rad=1, nd=1)

bvp_solver = beluga.bvp_algorithm('spbvp')

guess_maker = beluga.guess_generator(
    'auto',
    start=[h_0/meter, v_0/meter*second, m_0/kilogram],    # Starting values for states in order
    costate_guess=-0.1,
    control_guess=[pi/3],
    time_integrate=0.1,
)

continuation_steps = beluga.init_continuation()

continuation_steps.add_step() \
    .num_cases(5) \
    .const('v_f', 0)

continuation_steps.add_step() \
    .num_cases(5) \
    .const('m_f', 3.6e3/kilogram)

continuation_steps.add_step() \
    .num_cases(10) \
    .const('rho0', 1.225*meter**3/kilogram)

continuation_steps.add_step() \
    .num_cases(10, spacing='log') \
    .const('m_f', tar_m_f/kilogram)

continuation_steps.add_step() \
    .num_cases(10, spacing='log') \
    .const('eps2', 0.001)

continuation_steps.add_step() \
    .num_cases(10, spacing='log') \
    .const('eps', 0.001)

beluga.add_logger(logging_level=logging.DEBUG)

sol_set = beluga.solve(
    ocp=ocp,
    method='indirect',
    bvp_algorithm=bvp_solver,
    steps=continuation_steps,
    guess_generator=guess_maker,
    autoscale=False,
)

# Plot Results
sol = sol_set[-1][-1]

rho0 = 1.225

sol.y[:,0] = sol.y[:,0]*meter
sol.y[:,1] = sol.y[:,1]*meter/second
sol.y[:,2] = sol.y[:,2]*kilogram
sol.t = sol.t*second

tf = sol.t[-1]

thrustfun = interpolate.interp1d(sol.t, (np.sin(sol.u[:,0])+1)/2)

if not test_run:
    vessel.auto_pilot.target_pitch_and_heading(90, 90)
    vessel.auto_pilot.engage()
    vessel.control.throttle = float(thrustfun(0))

    obt_frame = vessel.orbit.body.non_rotating_reference_frame

    ut = conn.add_stream(getattr, conn.space_center, 'ut')
    altitude = conn.add_stream(getattr, vessel.flight(), 'mean_altitude')
    speed = conn.add_stream(getattr, vessel.flight(obt_frame), 'vertical_speed')
    drag = conn.add_stream(getattr, vessel.flight(), 'drag')

print('T-MINUS')

while countdown != 0:
    print(countdown)
    countdown -= 1
    time.sleep(1)

print('blastoff')

if not test_run:
    ut0 = ut()
    time0 = time.time()
    vessel.control.activate_next_stage()

    from beluga.ivpsol import Trajectory
    gamma = Trajectory()

    gamma.t = np.array([ut()-ut0])
    gamma.y = np.array([[altitude(), speed(), vessel.mass]])
    gamma.u = np.array([[thrustfun(0)]])
    D = np.array([drag()])
    T = np.array([vessel.thrust])

    while (ut() - ut0) < tf/3:
        vessel.control.throttle = float(thrustfun(ut()-ut0))
        gamma.t = np.hstack((gamma.t, np.array([ut()-ut0])))
        gamma.y = np.vstack((gamma.y, np.array([[altitude(), speed(), vessel.mass]])))
        gamma.u = np.vstack((gamma.u, thrustfun(ut()-ut0)))
        D = np.vstack((D, np.array(drag())))
        T = np.vstack((T, np.array(vessel.thrust)))
        time.sleep(0.1)

print('mission complete boss')

plt.figure(1)
plt.plot(sol.t, sol.y[:, 0])
if not test_run:
    plt.plot(gamma.t, gamma.y[:, 0], linestyle='--', label='Measured')
plt.xlabel('Time [s]')
plt.ylabel('Altitude [nd]')
plt.title('Altitude Profile')
plt.grid('on')
plt.legend()
plt.show()

plt.figure(2)
plt.plot(sol.t, sol.y[:, 1])
if not test_run:
    plt.plot(gamma.t, gamma.y[:, 1], linestyle='--', label='Measured')
plt.xlabel('Time [s]')
plt.ylabel('Velocity [nd]')
plt.title('Velocity Profile')
plt.grid('on')
plt.legend()
plt.show()

plt.figure(3)
plt.plot(sol.t, sol.y[:, 2])
if not test_run:
    plt.plot(gamma.t, gamma.y[:, 2], linestyle='--', label='Measured')
plt.xlabel('Time [s]')
plt.ylabel('Mass [nd]')
plt.title('Mass Profile')
plt.grid('on')
plt.legend()
plt.show()

plt.figure(4)
plt.plot(sol.t, (thrust_max_sea*np.exp(-sol.y[:,0]/H) + thrust_max_vac*(1-np.exp(-sol.y[:,0]/H)))*(np.sin(sol.u[:,0])+1)/2, label='Thrust')
if not test_run:
    pass
    plt.plot(gamma.t, T[:,0], label='Measured Thrust', linestyle='--')

plt.plot(sol.t, 1/2 * d * sol.y[:,1]**2 * rho0 * np.exp(-sol.y[:,0]/H) * pi * (A/2)**2, label='Drag')
if not test_run:
    plt.plot(gamma.t, -D[:,0], label='Measured Drag', linestyle='--')
plt.plot([0, tf], [drag_max, drag_max], color='k', linestyle='--')
plt.xlabel('Time [s]')
plt.ylabel('Force [N]')
plt.title('Forces')
plt.grid('on')
plt.legend()
plt.show()
