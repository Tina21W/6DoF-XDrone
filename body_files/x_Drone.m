function [properties,aerodynamics,initialization] = x_Drone()

    %% --------------- Properties of the craft - XDrone in this case -------------------
    % New XDrone quantities, formatted with the same downstream-facing field
    % names used by the old x_Drone.m file.

    % Geometry
    properties.Xzylo.CoG_pos = 33.041e-3;     % craft CoG from leading edge [m]
    properties.Xzylo.b = 131.6e-3;            % span [m]
    properties.Xzylo.c = 115.21e-3;           % chord [m]
    properties.Xzylo.percentage_CoG = properties.Xzylo.CoG_pos/properties.Xzylo.c;
    properties.Xzylo.mass = 99.61e-3;         % craft-only mass [kg]
    properties.Xzylo.S = properties.Xzylo.b*properties.Xzylo.c; % reference area [m^2]
    properties.Xzylo.aspect_ratio = properties.Xzylo.b/properties.Xzylo.c;

    % Complete system values from the new file
    properties.mass_total = 105.538e-3;       % total mass [kg]
    properties.CoG_pos_total = 31.648e-3;     % total CoG from leading edge [m]
    properties.percentage_CoG_total = properties.CoG_pos_total/properties.Xzylo.c;
    properties.Area = properties.Xzylo.S;
    properties.chord = properties.Xzylo.c;
    properties.span = properties.Xzylo.b;
    properties.aspect_ratio = properties.Xzylo.aspect_ratio;

    % Inertia quantities from the new file.  I_CoG is about each component's
    % own CoG; I is shifted to the total-system CoG using the parallel-axis theorem.
    Xzylo_CoG_vec = [properties.Xzylo.CoG_pos; 0; 0];
    total_CoG_vec = [properties.CoG_pos_total; 0; 0];
    properties.Xzylo.I_CoG = 1e-4 * [3.445  0       0
                                     0      2.0785  0
                                     0      0       2.0785];
    d_Xzylo = Xzylo_CoG_vec - total_CoG_vec;
    properties.Xzylo.I = properties.Xzylo.I_CoG + ...
        properties.Xzylo.mass*((d_Xzylo'*d_Xzylo)*eye(3) - (d_Xzylo*d_Xzylo'));

    %% ------------------------- Complete System Geometry -----------------------------
    properties.hub_radius = 18e-3;            % [m]
    properties.hub_height = 35.6e-3;          % [m]
    properties.no_blades = 3;
    properties.design_RPM = 900;             % [RPM]
    properties.design_speed = 15;             % [m/s]
    properties.f_coeff = 1.97e-2;              % friction torque coefficient

    %% ------------------------- Environmental Properties -----------------------------
    properties.rho = 1.225;
    properties.V_wind_i = [0; 0; 0];          % wind velocity in inertial frame [m/s]
    properties.g = 9.81;

    %% ------------------------- Propeller / Motor Group ------------------------------
    % New file gives the prop+motor assembly as one combined group.  Therefore
    % properties.Prop contains the combined assembly.  properties.Motor is kept
    % below as a zero-mass compatibility placeholder so old code that accesses
    % properties.Motor.* does not fail or double-count mass/inertia.

    % Negative RPM = counter-rotating relative to XDrone spin
    Prop_rpm = -12759;                        % [RPM]
    properties.Prop.omega = 2*pi/60 * Prop_rpm; % [rad/s]
    properties.Prop.mass = 5.928e-3;          % combined prop+motor mass [kg]
    properties.Prop.CoG_pos = [8.231e-3; 0; 0]; % combined prop+motor CoG [m]
    properties.Prop.radius = 3*25.4e-3/2;     % 3 inch diameter prop radius [m]

    properties.Prop.I_CoG = 1e-9 * [753  0    0
                                    0    500  0
                                    0    0    500]; % about prop+motor group CoG
    d_Prop = properties.Prop.CoG_pos - total_CoG_vec;
    properties.Prop.I = properties.Prop.I_CoG + ...
        properties.Prop.mass*((d_Prop'*d_Prop)*eye(3) - (d_Prop*d_Prop'));
    properties.Prop.I_yyzz = properties.Prop.I(2,2);

    % Propeller aerodynamic coefficients from the new file
    properties.Prop.k_T_0 = 0.205;
    properties.Prop.k_P_0 = properties.Prop.k_T_0/1.5;
    properties.Prop.k_Q_0 = properties.Prop.k_P_0/(2*pi);

    % Propeller speed controller -- restored from old file for compatibility
    properties.Prop.control.omega_Kp = 500;  % proportional gain [(rad/s)/(m/s)]
    properties.Prop.control.omega_Ki = 200;   % integral gain [(rad/s)/(m)]
    properties.Prop.control.int_min = -2;     % minimum integrated speed error [m]
    properties.Prop.control.int_max = 2;      % maximum integrated speed error [m]
    properties.Prop.control.omega_min = 0;    % minimum relative prop speed magnitude [rad/s]
    properties.Prop.control.omega_max = 4000; % maximum relative prop speed magnitude [rad/s]

    % Velocity direction controller -- restored from old file for compatibility
    properties.Direction_control.T_max = 0.05;      % maximum transverse control moment [Nm]
    properties.Direction_control.error_lim = pi/6; % direction error for maximum moment [rad]
    properties.Direction_control.Kd = 0.005;        % angular-rate damping gain [Nm/(rad/s)]

    % Motor compatibility block.  The new file already includes motor mass and
    % inertia in properties.Prop, so these fields are intentionally zeroed to
    % avoid double-counting if any later code sums component masses/inertias.
    properties.Motor.mass = 0;                % [kg]
    properties.Motor.radius = 0;              % [m]
    properties.Motor.L = 0;                   % [m]
    properties.Motor.CoG_pos = properties.Prop.CoG_pos;
    properties.Motor.I = zeros(3);


    %% ------------------------- Compatibility aliases -------------------------------
    % Minimal aliases for project code that builds sim from properties.Xzylo instead
    % of passing the full properties struct. These do not change the physical
    % quantities above; they only expose the same Prop/Motor/control fields in the
    % nested craft struct so calls such as sim.Prop.I keep working.
    properties.I = properties.Xzylo.I + properties.Prop.I + properties.Motor.I;
    properties.Xzylo.mass_total = properties.mass_total;
    properties.Xzylo.CoG_pos_total = properties.CoG_pos_total;
    properties.Xzylo.percentage_CoG_total = properties.percentage_CoG_total;
    properties.Xzylo.Area = properties.Area;
    properties.Xzylo.chord = properties.chord;
    properties.Xzylo.span = properties.span;
    properties.Xzylo.I_total = properties.I;
    properties.Xzylo.Prop = properties.Prop;
    properties.Xzylo.Motor = properties.Motor;
    properties.Xzylo.Direction_control = properties.Direction_control;

    %% Aerodynamics Group
    %----------------------- Aerodynamic coefficient functions Drone -------------------
    aerodynamics_xzylo                       % run aerodynamic file for interpolation of data
    aerodynamics.Xzylo.C_L = @(angle) C_L_interp(angle);
    aerodynamics.Xzylo.C_D = @(angle) C_D_interp(angle);
    aerodynamics.Xzylo.CoP_frac = @(angle) CoP_interp(angle);
    % aerodynamics.Xzylo.C_m = @(angle) (x_CoG - x_CoP(angle)) * ...
    %     (C_L*sin(angle) + C_D*cos(angle));
    % aerodynamics.Xzylo.C_Y = @(angle) 0;

    % Trim conditions
    [properties.alpha_trim, properties.V_trim, properties.Thrust_req] = calculate_trim(...
        properties.percentage_CoG_total, properties.mass_total*properties.g, properties.rho, ...
        properties.Area, aerodynamics.Xzylo.C_L, aerodynamics.Xzylo.C_D, aerodynamics.Xzylo.CoP_frac);
    properties.Prop.control.u_des = properties.V_trim; % desired body speed from trim [m/s]

    % ------------------------- Motor Blades Aerodynamics -----------------------------
    [V_GRID, RPM_GRID, Torque_Map, Drag_Map] = generate_turbine_maps(properties);
    aerodynamics.mblades.V_grid = V_GRID;
    aerodynamics.mblades.RPM_grid = RPM_GRID;
    aerodynamics.mblades.Torque_map = Torque_Map;
    aerodynamics.mblades.Drag_map = Drag_Map;

    % Optional external force/moment (written in BODY frame)
    aerodynamics.Fext = @(t) (t >= 0 && t <= 90) * [0.229; 0; 0];
    aerodynamics.Mext = @(t) (t >= 0 && t <= 30) * [0; 0; 0];

    % ------------------------- Initialization states ---------------------------------
    initialization.launch_angle = rad2deg(properties.alpha_trim);          % launch angle in degrees
    initialization.V_mag = properties.V_trim;               % magnitude of the launch velocity [m/s]
    initialization.Omega_mag = 12.85;            % rotational speed at launch [RPS]
    initialization.AoA = rad2deg(properties.alpha_trim);

    initialization.v0 = [initialization.V_mag*cosd(initialization.AoA), 0, ...
                         initialization.V_mag*sind(initialization.AoA)]; % [u v w] [m/s]
    initialization.omega0 = [2*pi*initialization.Omega_mag, 0, 0];        % [p q r] [rad/s]
    initialization.euler0 = [deg2rad(0), deg2rad(initialization.launch_angle), deg2rad(0)]; % [yaw pitch roll]
    initialization.quat0 = eul2quat(initialization.euler0);               % [w x y z]
    initialization.pos0 = [0 0 -600];                                      % inertial position [m]
    initialization.u_error_int0 = 0;                                      % integral speed error [m]
    initialization.tf = 20;                                               % maximum simulation time [s]

end
