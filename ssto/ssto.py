import time
import beluga
import krpc
import logging
from math import pi
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

test_run = False
countdown = 3

if not test_run:
    conn = krpc.connect(name='ksp-astrogator')
    vessel = conn.space_center.active_vessel
    print('Connected with ' + vessel.name)


ocp = beluga.OCP('ssto')

# Define independent variables
ocp.independent('t', 's')

# Define equations of motion
ocp.state('x', 'v_x', 'm') \
   .state('y', 'v_y', 'm') \
   .state('v_x', '(F_sea*exp(-y/Hscale) + F_vac*(1-exp(-y/Hscale)))/mass*cos(theta) - D/mass*v_x/sqrt(v_x**2 + v_y**2)', 'm/s') \
   .state('v_y', '(F_sea*exp(-y/Hscale) + F_vac*(1-exp(-y/Hscale)))/mass*sin(theta) - D/mass*v_y/sqrt(v_x**2 + v_y**2) - g', 'm/s') \
   .state('mass', 'md*eps', 'kg')

initial_mass = 10.665e3
Hscale = 5600
F_sea = 205.16e3
F_vac = 240e3

# Define controls
ocp.control('theta', 'rad')

ocp.quantity('D', '1/2*rho_ref*exp(-y/Hscale)*CD*A*(v_x**2 + v_y**2)')

# Define constants
ocp.constant('F_sea', 2.1e6, 'newton')
ocp.constant('F_vac', 2.1e6, 'newton')
ocp.constant('A', 7.069, 'm^2')
ocp.constant('mu', 3.986004e14, 'm^3/s^2')
ocp.constant('Re', 6378100, 'm')
ocp.constant('CD', 0.2, '1')
ocp.constant('rho_ref', 0, 'kg/m^3')
ocp.constant('Hscale', 8.44e3, 'm')
ocp.constant('g', 9.80665, 'm/s^2')
ocp.constant('md', -78.95, 'kg/s')
ocp.constant('eps', 0.000, '1')
ocp.constant('eps2', 10, 's')
ocp.constant('mupper', 65000, 'kg')
ocp.constant('mlower', 2000, 'kg')

ocp.constant('x_0', 0, 'm')
ocp.constant('y_0', 0, 'm')
ocp.constant('v_x_0', 0, 'm/s')
ocp.constant('v_y_0', 0.01, 'm/s')
ocp.constant('mass_0', 60880, 'kg')

ocp.constant('y_f', 1.8e5, 'm')
ocp.constant('v_y_f', 0, 'm/s')

# Define costs
ocp.path_cost('1', '1')
# ocp.terminal_cost('-mass', 'kg')

# Define constraints
ocp.constraints() \
    .initial('x - x_0', 'm')    \
    .initial('y - y_0', 'm') \
    .initial('v_x - v_x_0', 'm/s')  \
    .initial('v_y - v_y_0', 'm/s')  \
    .initial('mass - mass_0', 'kg') \
    .terminal('y - y_f', 'm') \
    .terminal('v_x - sqrt(mu/(y_f+Re))', 'm/s') \
    .terminal('v_y - v_y_f', 'm/s')

ocp.scale(m='y', s='y/v_x', kg='mass', newton='mass*v_x^2/y', rad=1)

bvp_solver = beluga.bvp_algorithm('spbvp')

guess_maker = beluga.guess_generator('auto',
                start=[0, 0, 0, 0.01, 60880],          # Starting values for states in order
                costate_guess = -0.1,
                control_guess=[0],
                use_control_guess=False
)

beluga.add_logger(logging_level=logging.DEBUG)

sol_set = beluga.solve(ocp=ocp,
             method='indirect',
             bvp_algorithm=bvp_solver,
             steps=None,
             guess_generator=guess_maker, autoscale=True)

guess_maker = beluga.guess_generator('static', solinit=sol_set[-1][-1])

continuation_steps = beluga.init_continuation()

continuation_steps.add_step('bisection') \
                .num_cases(10) \
                .const('F_sea', 205.16e3) \
                .const('F_vac', 240e3) \
                .const('Re', 600000) \
                .const('A', pi*(0.7**2)) \
                .const('mu', 3.5e12) \
                .const('mass_0', initial_mass) \
                .const('Hscale', 5600) \
                .const('y_f', 62.5e3)
# 240.16e3
continuation_steps.add_step('bisection') \
                .num_cases(10) \
                .const('eps', 1)

# continuation_steps.add_step('bisection') \
#                 .num_cases(10) \
#                 .const('mupper', initial_mass+3000) \
#                 .const('mlower', initial_mass-8000)

# continuation_steps.add_step('bisection') \
#                 .num_cases(10) \
#                 .const('rho_ref', 1.225)

sol_set = beluga.solve(ocp=ocp,
             method='indirect',
             bvp_algorithm=bvp_solver,
             steps=continuation_steps,
             guess_generator=guess_maker, autoscale=True)

sol = sol_set[-1][-1]

rho0 = 1.225

sol.y[:,0] = sol.y[:,0]
sol.y[:,1] = sol.y[:,1]
sol.y[:,2] = sol.y[:,2]
sol.y[:,3] = sol.y[:,3]
sol.y[:,4] = sol.y[:,4]
sol.t = sol.t

tf = sol.t[-1]

anglefun = interpolate.interp1d(sol.t, sol.u[:,0]*180/np.pi)

if not test_run:
    vessel.auto_pilot.target_pitch_and_heading(float(anglefun(0)), 90)
    vessel.auto_pilot.engage()
    vessel.control.throttle = 1.0

    obt_frame = vessel.orbit.body.non_rotating_reference_frame

    ut = conn.add_stream(getattr, conn.space_center, 'ut')
    altitude = conn.add_stream(getattr, vessel.flight(), 'mean_altitude')
    speed = conn.add_stream(getattr, vessel.flight(obt_frame), 'vertical_speed')
    drag = conn.add_stream(getattr, vessel.flight(), 'drag')
    pitch = conn.add_stream(getattr, vessel.flight(), 'pitch')

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
    gamma.u = np.array([[anglefun(0)]])
    D = np.array([drag()])
    T = np.array([vessel.thrust])
    P = np.array([pitch()])

    """
    Main Running Loop
    """
    while (ut() - ut0) < tf:
        vessel.auto_pilot.target_pitch_and_heading(float(anglefun(ut()-ut0)), 90)
        gamma.t = np.hstack((gamma.t, np.array([ut()-ut0])))
        gamma.y = np.vstack((gamma.y, np.array([[altitude(), speed(), vessel.mass]])))
        gamma.u = np.vstack((gamma.u, anglefun(ut()-ut0)))
        D = np.vstack((D, np.array(drag())))
        T = np.vstack((T, np.array(vessel.thrust)))
        P = np.hstack((P, pitch()))
        time.sleep(0.1)

    vessel.control.throttle = 0.0

print('mission complete boss')

plt.figure(1)
plt.plot(sol.t, sol.y[:, 1])
if not test_run:
    plt.plot(gamma.t, gamma.y[:, 0], linestyle='--', label='Measured')
plt.xlabel('Time [s]')
plt.ylabel('Altitude [m]')
plt.title('Altitude Profile')
plt.grid('on')
plt.legend()
plt.show()

plt.figure(2)
plt.plot(sol.t, sol.y[:, 3])
if not test_run:
    plt.plot(gamma.t, gamma.y[:, 1], linestyle='--', label='Measured')
plt.xlabel('Time [s]')
plt.ylabel('Vertical Velocity [m/s]')
plt.title('Velocity Profile')
plt.grid('on')
plt.legend()
plt.show()

plt.figure(3)
plt.plot(sol.t, sol.y[:, 4])
if not test_run:
    plt.plot(gamma.t, gamma.y[:, 2], linestyle='--', label='Measured')
plt.xlabel('Time [s]')
plt.ylabel('Mass [kg]')
plt.title('Mass Profile')
plt.grid('on')
plt.legend()
plt.show()

plt.figure(4)
plt.plot(sol.t, sol.u[:,0]*180/np.pi, label='Pitch')
if not test_run:
    plt.plot(gamma.t, P, label='Measured Pitch', linestyle='--')
plt.xlabel('Time [s]')
plt.ylabel('Pitch [deg]')
plt.title('Steering Angle')
plt.grid('on')
plt.legend()
plt.show()

plt.plot()
plt.plot(sol.t, F_sea*np.exp(-sol.y[:,1]/Hscale) + F_vac*(1-np.exp(-sol.y[:,1]/Hscale)), label='Thrust')
if not test_run:
    plt.plot(gamma.t, T[:,0], label='Measured Thrust', linestyle='--')
plt.xlabel('Time [s]')
plt.ylabel('Thrust [N]')
plt.title('Thrust')
plt.grid('on')
plt.legend()
plt.show()
