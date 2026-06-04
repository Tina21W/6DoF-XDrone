function [target_waypoint, current_wp_idx] = waypoint_manager(pos_i, ctrl)
    % Manages sequential waypoint tracking and returns next waypoint coordinate
    % Inputs: pos_i = inertial position [3x1]
    %         ctrl = control struct with waypoints, waypoint_radius, loop options
    % Outputs:
    %   target_waypoint = coordinate of next waypoint [3x1]
    %   current_wp_idx = index of current waypoint
    
    persistent wp_idx last_wp_idx print_counter;
    
    % Initialize persistent variables
    if isempty(wp_idx)
        wp_idx = 1;
        last_wp_idx = 0;
        print_counter = 0;
    end
    
    % Extract and normalize waypoints
    if isfield(ctrl, 'waypoints')
        waypoints = ctrl.waypoints;
    else
        waypoints = ctrl.waypoint;
    end
    
    if size(waypoints, 1) ~= 3
        waypoints = waypoints';
    end
    
    N = size(waypoints, 2);
    
    % Get waypoint radius (default 10m)
    radius = 10;
    if isfield(ctrl, 'waypoint_radius')
        radius = ctrl.waypoint_radius;
    end
    
    % Check if we've reached current waypoint
    dist_to_current = norm(waypoints(:, wp_idx) - pos_i);
    if dist_to_current < radius
        if isfield(ctrl, 'loop') && ctrl.loop
            wp_idx = mod(wp_idx, N) + 1;
        else
            wp_idx = min(wp_idx + 1, N);
        end
    end
    
    % Log waypoint switches
    if wp_idx ~= last_wp_idx
        fprintf('[WP] Switched to WP%d/%d | pos=[%.1f, %.1f, %.1f] | dist_prev=%.1f m\n', ...
            wp_idx, N, pos_i(1), pos_i(2), pos_i(3), dist_to_current);
        last_wp_idx = wp_idx;
    end
    
    % Periodic logging
    print_counter = print_counter + 1;
    if mod(print_counter, 200) == 0
        d = norm(waypoints(:, wp_idx) - pos_i);
        fprintf('[WP] Targeting WP%d/%d | pos=[%.1f, %.1f, %.1f] | dist=%.1f m\n', ...
            wp_idx, N, pos_i(1), pos_i(2), pos_i(3), d);
    end
    
    % Return waypoint coordinate
    target_waypoint = waypoints(:, wp_idx);
    current_wp_idx = wp_idx;
end
