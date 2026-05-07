function [T_vector] = controller(pos_i, v_b, waypoint, R_ib)
    % 3D Waypoint tracking control
    % Inputs: pos_i = inertial position, v_b = body-frame velocity
    %         waypoint = [x_wp; y_wp; z_wp], R_ib = inertial->body rotation
    % Output: T_vector = [0; My; Mz] - torque in body frame y-z plane

    wp = waypoint;

    % Desired direction in inertial frame
    desired_dir = wp - pos_i;
    desired_dir = desired_dir / norm(desired_dir);

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
    T_max = 0.015;         % tune this
    error_lim = pi/2;      % tune this

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

    T_vector = [0; Ty; Tz];

    %fprintf('pos_i = [%g, %g, %g]\n', pos_i(1), pos_i(2), pos_i(3));
    %fprintf('wp = [%g, %g, %g]\n', wp(1), wp(2), wp(3));
    %fprintf('v_b = [%g, %g, %g]\n', v_b(1), v_b(2), v_b(3));
    %fprintf('v_i = [%g, %g, %g]\n', v_i(1), v_i(2), v_i(3));
    %fprintf('desired_dir = [%g, %g, %g]\n', desired_dir(1), desired_dir(2), desired_dir(3));
    %fprintf('current_dir = [%g, %g, %g]\n', current_dir(1), current_dir(2), current_dir(3));
    %fprintf('delta = %g\n', delta);

end
