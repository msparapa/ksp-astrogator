import beluga
import logging
import matplotlib.pyplot as plt
import numpy as np
import krpc
import time
from scipy import interpolate

test_run = False
countdown = 3

"""
Kerbin params
"""

Hscale = 5600
mu = 3.5316e12
Re = 600000

twostage_thrust0_sea = 374.19e3
twostage_thrust0_vac = 400e3
twostage_thrust1_sea = 205.16e3
twostage_thrust1_vac = 240e3

twostage_mass0 = 29.365e3
twostage_mass0f = twostage_mass0 - 18.8e3
twostage_mass1 = twostage_mass0 - 2*10.125e3 - 2e3 - 0.9e3

twostage_massflow0 = -131.58
twostage_massflow1 = -78.95

y_f = 100e3


if not test_run:
    conn = krpc.connect(name='ksp-astrogator')
    vessel = conn.space_center.active_vessel
    print('Connected with ' + vessel.name)

ocp = beluga.OCP('twostage')

# Define independent variables
ocp.independent('t', 's')

# Define equations of motion
ocp.state('x', 'v_x', 'm') \
   .state('y', 'v_y', 'm') \
   .state('v_x', 'Thrust/current_mass*cos(theta) - D/current_mass*v_x/sqrt(v_x**2 + v_y**2)', 'm/s') \
   .state('v_y', 'Thrust/current_mass*sin(theta) - D/current_mass*v_y/sqrt(v_x**2 + v_y**2) - g', 'm/s') \
   .state('mass', 'mass_flow*eps', 'kg')

# Define controls
ocp.control('theta', 'rad')

ocp.quantity('engine0', 'F0_sea*exp(-y/Hscale) + F0_vac*(1-exp(-y/Hscale))')
ocp.quantity('engine1', 'F1_sea*exp(-y/Hscale) + F1_vac*(1-exp(-y/Hscale))')

ocp.quantity('Thrust', 'engine0/(1+exp(rash*(mass_0f - mass))) + engine1/(1+exp(rash*(mass - mass_0f)))')
ocp.quantity('mass_flow', 'md0/(1+exp(rash*(mass_0f - mass))) + md1/(1+exp(rash*(mass - mass_0f)))')
ocp.quantity('current_mass', 'mass - drop_mass/(1+exp(rash*(mass - mass_0f)))')
ocp.quantity('D', '1/2*rho_ref*exp(-y/Hscale)*CD*A*(v_x**2 + v_y**2)')


# Define constants
ocp.constant('F0_sea', 2.1e6, 'newton')
ocp.constant('F0_vac', 2.1e6, 'newton')
ocp.constant('F1_sea', 2.1e6, 'newton')
ocp.constant('F1_vac', 2.1e6, 'newton')
ocp.constant('A', np.pi*(1.25/2)**2, 'm^2')
ocp.constant('mu', 3.986004e14, 'm^3/s^2')
ocp.constant('Re', 6378100, 'm')
ocp.constant('CD', 0.1, '1')
ocp.constant('rho_ref', 0, 'kg/m^3')
ocp.constant('Hscale', 8.44e3, 'm')
ocp.constant('g', 9.80665, 'm/s^2')
ocp.constant('md0', -807.6, 'kg/s')
ocp.constant('md1', -807.6, 'kg/s')
ocp.constant('eps', 0.000, '1')

ocp.constant('x_0', 0, 'm')
ocp.constant('y_0', 0, 'm')
ocp.constant('v_x_0', 0, 'm/s')
ocp.constant('v_y_0', 0.01, 'm/s')
ocp.constant('mass_0', 60880, 'kg')
ocp.constant('mass_0f', 0, 'kg')
ocp.constant('drop_mass', 0, 'kg')

ocp.constant('rash', 0.01, '1')

ocp.constant('y_f', 1.5e5, 'm')
ocp.constant('v_y_f', 0, 'm/s')

# Define costs
ocp.path_cost('1', '1')

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
                .num_cases(20) \
                .const('eps', 1) \
                .const('mu', mu) \
                .const('y_f', y_f) \
                .const('Hscale', Hscale) \
                .const('Re', Re) \
                .const('mass_0', twostage_mass0) \
                .const('mass_0f', twostage_mass0f) \
                .const('F0_sea', twostage_thrust0_sea) \
                .const('F0_vac', twostage_thrust0_vac) \
                .const('F1_sea', twostage_thrust1_sea) \
                .const('F1_vac', twostage_thrust1_vac) \
                .const('md0', twostage_massflow0) \
                .const('md1', twostage_massflow1) \
                .const('drop_mass', twostage_mass0f - twostage_mass1) \
                .const('rash', 1)

continuation_steps.add_step('bisection') \
                .num_cases(10) \
                .const('rho_ref', 1.225)

continuation_steps.add_step('bisection') \
                .num_cases(10) \
                .const('rash', 100)

try:
    sol_set = beluga.solve(ocp=ocp,
                 method='indirect',
                 bvp_algorithm=bvp_solver,
                 steps=continuation_steps,
                 guess_generator=guess_maker, autoscale=True)
except:
    print('snake? SNAAAAAAKE!!!!!')

sol = sol_set[-1][-1]
tf = sol.t[-1]

massfun = interpolate.interp1d(sol.t, sol.y[:,4])
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
    stage = 1

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
        if massfun(ut()-ut0) <= twostage_mass0f and stage == 1:
            print('next stage')
            vessel.control.activate_next_stage()
            time.sleep(0.1)
            vessel.control.activate_next_stage()
            stage += 1

        time.sleep(0.1)

    vessel.control.throttle = 0.0

print('mission complete boss')


plt.plot(sol.y[:,0]/1000, sol.y[:,1]/1000)
plt.xlabel('Downrange [km]')
plt.ylabel('Altitude [km]')
plt.title('Time Optimal Launch Trajectory')
plt.grid('on')
plt.show()

plt.plot(sol.t, sol.y[:,2], label='Horizontal')
plt.plot(sol.t, sol.y[:,3], label='Vertical')
if not test_run:
    plt.plot(gamma.t, gamma.y[:, 1], linestyle='--', label='Measured Vertical')
plt.xlabel('Time [s]')
plt.ylabel('Velocity [km/s]')
plt.title('Velocities')
plt.legend()
plt.grid('on')
plt.show()

plt.plot(sol.t, sol.u*180/np.pi)
if not test_run:
    plt.plot(gamma.t, P, label='Measured Pitch', linestyle='--')
plt.xlabel('Time [s]')
plt.ylabel('Control [degrees]')
plt.title('Steering Angle')
plt.grid('on')
plt.show()

plt.plot(sol.t, sol.y[:,4] - (twostage_mass0f - twostage_mass1)/(1+np.exp(sol.aux['const']['rash']*(sol.y[:,4] - twostage_mass0f))))
if not test_run:
    plt.plot(gamma.t, gamma.y[:, 2], linestyle='--', label='Measured')
plt.plot([sol.t[0], sol.t[-1]], [twostage_mass0, twostage_mass0], linestyle='--', color='k')
plt.plot([sol.t[0], sol.t[-1]], [twostage_mass0f, twostage_mass0f], linestyle='--', color='k')
plt.xlabel('Time [s]')
plt.ylabel('Mass [kg]')
plt.title('2 Stage Mass')
plt.grid('on')
plt.show()

engine0 = twostage_thrust0_sea*np.exp(-sol.y[:,1]/Hscale) + twostage_thrust0_vac*(1-np.exp(-sol.y[:,1]/Hscale))
engine1 = twostage_thrust1_sea*np.exp(-sol.y[:,1]/Hscale) + twostage_thrust1_vac*(1-np.exp(-sol.y[:,1]/Hscale))
thrust = engine0/(1+np.exp(sol.aux['const']['rash']*(twostage_mass0f - sol.y[:,4]))) + engine1/(1+np.exp(sol.aux['const']['rash']*(sol.y[:,4] - twostage_mass0f)))
plt.plot(sol.t, thrust)
if not test_run:
    plt.plot(gamma.t, T[:, 0], linestyle='--', label='Measured')
plt.xlabel('Time [s]')
plt.ylabel('Thrust [newtons]')
plt.title('Thrust')
plt.grid('on')
plt.show()
