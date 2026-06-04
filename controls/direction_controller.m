function [T_mag, T_mag_undamped, alpha, beta, delta] = direction_controller(pos_i, v_b, omega_b, ctrl, R_ib, sim, q); 
    [target_waypoint, current_wp_idx] = waypoint_manager(pos_i, ctrl);
    % Pure direction-tracking PD controller (no waypoint management)
    % Inputs: pos_i = inertial position [3x1]
    %         target_waypoint = waypoint coordinate in inertial frame [3x1]
    %         v_b = body-frame velocity [3x1]
    %         omega_b = body angular rate [3x1]
    %         R_ib = inertial->body rotation matrix [3x3]
    %         sim = simulation struct with Direction_control parameters
    % Outputs:
    %   T_cmd = [0; My; Mz] - torque command in body frame [3x1]
    %   T_mag = magnitude of torque (after damping)
    %   alpha = desired direction angle in de-spun yz-plane
    %   beta = direction error angle in de-spun yz-plane
    %   delta = angular heading error (desired vs current direction)

    % Compute desired direction from position to target
    to_target = target_waypoint - pos_i;
    dist_to_target = norm(to_target);
    
    if dist_to_target < 1e-6
        % At waypoint, no control needed
        desired_dir = [1; 0; 0];
    else
        desired_dir = to_target / dist_to_target;
    end

    % Current direction - convert body velocity to inertial
    v_i = R_ib' * v_b;
    vel_mag = norm(v_i);
    if vel_mag < 1e-10
        current_dir = [1; 0; 0];
    else
        current_dir = v_i / vel_mag;
    end

    % Compute desired direction in body frame
    desired_b = R_ib * desired_dir;

    % Rotate desired body vector into the de-spun frame
    %eul = quat2eul(q(:)', 'ZYX');
    %phi = eul(3);
    % Extract roll angle (phi) directly from R_ib matrix
    phi = atan2(-R_ib(2,3), R_ib(3,3));
    c_phi = cos(phi);
    s_phi = sin(phi);
    R_x_phi = [1 0 0; 0 c_phi s_phi; 0 -s_phi c_phi];
    desired_ds = R_x_phi * desired_b;
    
    % Angular error (angle between desired and current)
    delta = acos(max(-1, min(1, dot(desired_dir, current_dir))));
    
    % Vector error
    error_i = desired_dir - current_dir;
    
    % Transform error to body frame
    error_b = R_ib * error_i;

    % Transform error to despun body frame
    error_ds = R_x_phi * error_b;

    % Angle in yz-plane relative to de-spun z-axis for desired direction
    alpha = atan2(desired_ds(2), desired_ds(3));

    % Angle in yz-plane relative to de-spun z-axis for vector error
    beta = atan2(error_ds(2), error_ds(3));


    % Extract control parameters
    if isfield(sim.prop, 'Direction_control')
        T_max = sim.prop.Direction_control.T_max;
        error_lim = sim.prop.Direction_control.error_lim;
        Kd = sim.prop.Direction_control.Kd;
    else
        T_max = 0.05;      % maximum torque
        error_lim = pi/2;
        Kd = 0.002;
    end

    G = T_max/error_lim;   % Proportional gain

    % Compute proportional command based on heading error
    if delta > error_lim
        T_mag = T_max;
    else
        T_mag = G * delta;
    end

    T_mag_undamped = T_mag;

    % Calculate perpendicular angular velocity magnitude for damping
    omega_perp = sqrt(omega_b(2)^2 + omega_b(3)^2);
    T_mag_damped = T_mag - Kd * omega_perp;
    T_mag_damped = max(0, T_mag_damped);

    T_mag = min(T_mag_damped, T_max);  % Return damped magnitude

end
