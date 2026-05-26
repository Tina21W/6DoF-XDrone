function [T_cmd, T_mag] = controller(pos_i, v_b, omega_b, ctrl, R_ib, sim)
    % 3D velocity-direction PD control for sequential waypoint tracking
    % Inputs: pos_i = inertial position, v_b = body-frame velocity
    %         omega_b = body angular rate, ctrl = control options struct,
    %         R_ib = inertial->body rotation, sim = simulation struct
    % Outputs:
    %   T_cmd = [0; My; Mz] - torque in body frame y-z plane
    %   T_mag    = magnitude of the torque vector before any OAP scaling

    persistent wp_idx last_wp_idx print_counter;

    if isfield(ctrl, 'waypoints')
        waypoints = ctrl.waypoints;
    else
        waypoints = ctrl.waypoint;
    end

    if size(waypoints, 1) ~= 3
        waypoints = waypoints';
    end

    N = size(waypoints, 2);
    radius = 10;
    if isfield(ctrl, 'waypoint_radius')
        radius = ctrl.waypoint_radius;
    end

    if isempty(wp_idx) || wp_idx < 1 || wp_idx > N
        wp_idx = 1;
        last_wp_idx = 0;
        print_counter = 0;
    end

    dist_to_current = norm(waypoints(:, wp_idx) - pos_i);
    if dist_to_current < radius
        if isfield(ctrl, 'loop') && ctrl.loop
            wp_idx = mod(wp_idx, N) + 1;
        else
            wp_idx = min(wp_idx + 1, N);
        end
    end

    if wp_idx ~= last_wp_idx
        fprintf('[WP] Switched to WP%d/%d | pos=[%.1f, %.1f, %.1f] | dist_prev=%.1f m\n', ...
            wp_idx, N, pos_i(1), pos_i(2), pos_i(3), dist_to_current);
        last_wp_idx = wp_idx;
    end

    print_counter = print_counter + 1;
    if mod(print_counter, 200) == 0
        d = norm(waypoints(:, wp_idx) - pos_i);
        fprintf('[WP] Targeting WP%d/%d | pos=[%.1f, %.1f, %.1f] | dist=%.1f m\n', ...
            wp_idx, N, pos_i(1), pos_i(2), pos_i(3), d);
    end

    wp = waypoints(:, wp_idx);

    % Desired direction in inertial frame
    to_wp = wp - pos_i;
    dist = norm(to_wp);
    if dist < 1e-6
        T_cmd = [0; 0; 0];
        return;
    end
    desired_dir = to_wp / dist;

    % Current direction - convert body velocity to inertial
    v_i = R_ib' * v_b;
    vel_mag = norm(v_i);
    if vel_mag < 1e-10
        current_dir = [1; 0; 0];
    else
        current_dir = v_i / vel_mag;
    end
    
    % Angular error (angle between desired and current)
    delta = acos(max(-1, min(1, dot(desired_dir, current_dir))));

    % Error as vector 
    error_i = desired_dir - current_dir;
    
    % Transform error to body frame
    error_b = R_ib * error_i;

    % Control parameters
    if isfield(sim.prop, 'Direction_control')
        T_max = sim.prop.Direction_control.T_max;
        error_lim = sim.prop.Direction_control.error_lim;
        Kd = sim.prop.Direction_control.Kd;
    else
        T_max = 0.05;      % maximum torque
        error_lim = pi/2;
        Kd = 0.002;
    end

    G = T_max/error_lim;   % Gain

    if delta > error_lim
        T_mag = T_max;
    else
        T_mag = G * delta;
    end

    if norm(error_b(2:3)) > 1e-10
        Ty = T_mag * error_b(2) / norm(error_b(2:3));
        Tz = T_mag * error_b(3) / norm(error_b(2:3));
    else
        Ty = 0;
        Tz = 0;
    end

    % Angular-rate damping turns the proportional direction controller into PD.
    Ty = Ty - Kd * omega_b(2);
    Tz = Tz - Kd * omega_b(3);

    Ty = min(max(Ty, -T_max), T_max);
    Tz = min(max(Tz, -T_max), T_max);

    T_cmd = [0; Ty; Tz];
    T_cmdi = R_ib' * T_cmd;

    %fprintf('pos_i = [%g, %g, %g]\n', pos_i(1), pos_i(2), pos_i(3));
    %fprintf('wp = [%g, %g, %g]\n', wp(1), wp(2), wp(3));
    %fprintf('v_b = [%g, %g, %g]\n', v_b(1), v_b(2), v_b(3));
    %fprintf('v_i = [%g, %g, %g]\n', v_i(1), v_i(2), v_i(3));
    %fprintf('desired_dir = [%g, %g, %g]\n', desired_dir(1), desired_dir(2), desired_dir(3));
    %fprintf('current_dir = [%g, %g, %g]\n', current_dir(1), current_dir(2), current_dir(3));
    %fprintf('delta = %g\n', delta);

end
