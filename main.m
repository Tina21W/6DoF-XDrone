%% ===============================
%% MAIN.M - 6DoF SIMULATION
%% ===============================

close all; clear; clc;
% globalData % load the global struct for logging data from utils_calc folder

%% ------------- INITIALIZATION -------------------------------------------
% Select vehicle configuration
vehicle_config = @x_Drone_0;  % Put the name of the body_file after @

% Load configuration
% ============================================
[sim.prop, sim.aero, sim.initial] = vehicle_config(); %properties, aerodynamics, initial conditions
% =============================================

x0 = [sim.initial.v0 sim.initial.omega0 sim.initial.quat0 sim.initial.pos0];
sim.initial.x0 = x0;

% --- Simulation setup ---
tspan = [0, sim.initial.tf];  % total time

%% ------------- SIMULATION HEADER ----------------------------------------
fprintf('==================================\n')
fprintf('   ANNULAR WING 6DOF SIMULATION\n')
fprintf('==================================\n')
fprintf('Vehicle configuration: %s\n', func2str(vehicle_config))
fprintf('----------------------------------\n')
fprintf('Launch Angle: %.1f°\n', sim.initial.launch_angle)
fprintf('Launch AoA: %.1f°\n', sim.initial.AoA)
fprintf('Launch Velocity: %.1f m/s\n', sim.initial.V_mag)
fprintf('Launch Rotational Speed: %.1f Hz\n', sim.initial.Omega_mag)
fprintf('----------------------------------\n')
fprintf('Trim AoA: %.2f°\n', rad2deg(sim.prop.alpha_trim))
fprintf('Trim Velocity: %.2f m/s\n', sim.prop.V_trim)
fprintf('Trim Thrust: %.3f N\n', sim.prop.Thrust_req)
fprintf('----------------------------------\n')

%% ------------- CONTROL OPTIONS -----------------------------------------
sim.options.control.control_law = @OAP;    % choose @OAP or @SPL

%% ------------- TRAJECTORY DEFINITION -----------------------------------
% Figure-8 loop: 4 waypoints tracing a horizontal figure-8.
% The drone crosses the centre twice per lap (once per lobe).
% Scale the whole shape by changing R.
% 12-point figure-8: two tangent circles, correctly ordered.
%
% Drone launches at [0,0] going +x. Right circle centre is at [0,+R].
% The drone enters at the bottom of the circle going +x (tangent point),
% then curves counter-clockwise. WP1 is the first arc point it reaches,
% only 22.5° off the launch heading — no abrupt first turn.
%
% All turns <= 45°, all radii >= 50m (physical min ~26m at 20m/s).
R  = 50;   % circle radius [m]
z0 = -15;  % altitude — must match initialization.pos0(3)

sim.options.control.waypoints = [
% right left
    50,  50,  z0;
    50,  100,  z0;
    

% % straight
%       500,  0, z0;

]';

% Min gap 38m, max turn 45°, crossing segments 66m — within gyro limits.
sim.options.control.waypoint_radius = 10.0;

% Loop forever: wrap wp4 back to wp1
sim.options.control.loop = false;

%% ------------- PLOTTING OPTIONS -----------------------------------------
sim.options.propeller_on = false;
sim.options.mBlades_on = false;

%% ------------- SPEED-HOLD (virtual propeller) ---------------------------
% Applies just enough thrust at each timestep to maintain a constant speed.
% This simulates a propeller keeping the drone at trim speed indefinitely.
% Turn on gravity in original_xzylo.m (g = 9.81) and this will compensate.
sim.options.speed_hold = true;
sim.options.V_target   = 20;   % target speed [m/s] — match V_mag in body file

sim.options.body_plotting = true;
sim.options.non_rotating_plotting = true;
sim.options.test_plotting = false;
sim.options.live_plotting_calculation = false; 
sim.options.radius_visualization = 5;

% ========================= INTEGRATOR ===============================
options_integration = options_premaker(sim.options.live_plotting_calculation, sim);

clear controller;  % reset persistent waypoint index before each run
fprintf('Starting integration...\n'); tic; 
% Start the integration process using the ODE solver
[t, x] = ode45(@(t,x) sixDoF_wrapper(t,x,sim), tspan, x0, options_integration);

time_ode = toc; 
fprintf('Integration finished in %.3f seconds.\n', time_ode);
% ====================================================================

%% ------------- POST-PROCESSING ------------------------------------------
tic;
SIM_DATA = compute_logged_data(t, x, sim);
time_log = toc;
fprintf('Post-processing finished in %.3f seconds.\n', time_log);

% SIM_DATA is already in the workspace!
plot_results(t, x, sim, SIM_DATA);
plot_nutation_error(t, x, sim);
plot_animation(t, x, sim);

% visualization_XZylo_orientation(sim);