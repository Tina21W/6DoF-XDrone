%% ===============================
%% MAIN.M - 6DoF SIMULATION
%% ===============================

close all; clear; clc;
% globalData % load the global struct for logging data from utils_calc folder

%% ------------- INITIALIZATION -------------------------------------------
% Select vehicle configuration
vehicle_config = @x_Drone;  % Put the name of the body_file after @

% Load configuration
% ============================================
[sim.prop, sim.aero, sim.initial] = vehicle_config(); %properties, aerodynamics, initial conditions
% =============================================

x0 = [sim.initial.v0 sim.initial.omega0 sim.initial.quat0 sim.initial.pos0 sim.initial.u_error_int0];
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
% 12-point figure-8: two tangent circles, ordered for smooth entry.
R  = 100;                    % circle radius [m]
z0 = sim.initial.pos0(3);   % altitude [m]

sim.options.control.waypoints = [
    R*cosd(315),  R + R*sind(315),  z0
    R*cosd(  0),  R + R*sind(  0),  z0
    R*cosd( 45),  R + R*sind( 45),  z0
    R*cosd( 90),  R + R*sind( 90),  z0
    R*cosd(135),  R + R*sind(135),  z0
    R*cosd(180),  R + R*sind(180),  z0
    R*cosd(135), -R + R*sind(135),  z0
    R*cosd(180), -R + R*sind(180),  z0
    R*cosd(225), -R + R*sind(225),  z0
    R*cosd(270), -R + R*sind(270),  z0
    R*cosd(315), -R + R*sind(315),  z0
    R*cosd(  0), -R + R*sind(  0),  z0
]';

sim.options.control.waypoint_radius = 20.0;
sim.options.control.loop = true;

%% ------------- PLOTTING OPTIONS -----------------------------------------
sim.options.propeller_on = true;
sim.options.mBlades_on = true;

sim.options.body_plotting = true;
sim.options.non_rotating_plotting = true;
sim.options.test_plotting = false;
sim.options.live_plotting_calculation = false; 
sim.options.radius_visualization = 5;

% ========================= INTEGRATOR ===============================
options_integration = options_premaker(sim.options.live_plotting_calculation, sim);

fprintf('Starting integration...\n'); tic; 
% Start the integration process using the ODE solver
clear controller;
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
plot_animation(t, x, sim);

% visualization_XZylo_orientation(sim);
