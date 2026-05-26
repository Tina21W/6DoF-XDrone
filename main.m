%% ===============================
%% MAIN.M - 6DoF SIMULATION
%% ===============================

close all; clear; clc;
% globalData % load the global struct for logging data from utils_calc folder

%% ------------- INITIALIZATION -------------------------------------------
% Select vehicle configuration
vehicle_config = @x_Drone;  % Put the name of the body_file after

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
% >>>  CHANGE THIS ONE LINE to switch trajectory  <<<
%      Options: 'figure8' | 'hstep' | 'vstep' | 'straight'
trajectory = 'straight';

R  = 45;   % circle radius [m]  — scales figure8 and circle
z0 = -50;  % altitude [m, NED]  — must match sim.initial.pos0(3)

loop = true;  % whether to cycle through waypoints continuously

switch trajectory

    case 'figure8'
        % 12-point horizontal figure-8.  Two CCW tangent circles.
        % Drone launches at [0,0] heading +x; right lobe centre [0,+R].
        % 270° (right lobe) and 90° (left lobe) are the crossing points — skipped.
        wp = [
        %  Right lobe CCW from entry (centre [0,+R])
            R*cosd(315),  R + R*sind(315),  z0;   % wp1  — 22.5° off launch heading
            R*cosd(  0),  R + R*sind(  0),  z0;   % wp2  — right tip
            R*cosd( 45),  R + R*sind( 45),  z0;   % wp3
            R*cosd( 90),  R + R*sind( 90),  z0;   % wp4  — top
            R*cosd(135),  R + R*sind(135),  z0;   % wp5
            R*cosd(180),  R + R*sind(180),  z0;   % wp6
        %  Left lobe CCW from re-entry (centre [0,-R])
            R*cosd(  0), -R + R*sind(  0),  z0;   % wp7
            R*cosd(315), -R + R*sind(315),  z0;   % wp8
            R*cosd(270), -R + R*sind(270),  z0;   % wp9  — bottom
            R*cosd(225), -R + R*sind(225),  z0;   % wp10
            R*cosd(180), -R + R*sind(180),  z0;   % wp11 — left tip
            R*cosd(135), -R + R*sind(135),  z0;   % wp12
        ];

    case 'hstep'
        % Horizontal (lateral) step.
        % Fly +x 100 m, side-step +y 50 m, resume +x 100 m at new offset.
        % Tests lateral position response at constant altitude.
        wp = [
            100,   0,  z0;   % hs1 — end of initial straight leg
            100,  150,  z0;   % hs2 — end of side-step  (step input)
            500,  150,  z0;   % hs3 — resume forward at new lateral offset
            
        ];

    case 'vstep'
        % Vertical (altitude) step.
        % Fly +x 100 m, climb 15 m, resume +x 100 m at new altitude.
        % NED convention: more-negative z = higher altitude.
        wp = [
            100,   0,   z0;   % vs1 — end of initial leg   (15 m AGL)
            100,   0,  -250;   % vs2 — top of climb         (30 m AGL, +15 m step)
            350,   0,  -250;   % vs3 — cruise at new altitude
            350,   0,   z0;
            500,   0,   z0;
        ];

    case 'straight'
        % Straight-line trajectory at constant altitude.
        % Fly +x from the start position in a straight line.
        wp = [
            100,   0,  z0;   % end of first straight leg
            500,   0,  z0;   % end of second straight leg
            900,   0,  z0;   % extend the straight line further
        ];
        loop = false;

    otherwise
        error('Unknown trajectory ''%s''. Choose: figure8 | hstep | vstep | straight', trajectory);
end

sim.options.control.waypoints = wp';

% Min gap 38m, max turn 45°, crossing segments 66m — within gyro limits.
sim.options.control.waypoint_radius = 10.0;

% Loop mode: only true for repeating trajectories
sim.options.control.loop = loop;

%% ------------- PLOTTING OPTIONS -----------------------------------------
sim.options.propeller_on = true;
sim.options.mBlades_on = true;

sim.options.body_plotting = true;
sim.options.non_rotating_plotting = true;
sim.options.test_plotting = true;
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
