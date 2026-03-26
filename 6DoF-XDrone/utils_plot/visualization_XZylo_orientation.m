function visualization_XZylo_orientation(params)
    global SIM_DATA;

    % --- Setup ---
    force_scale = 5;
    h_fig = figure('Name', 'Cylinder Orientation with θ-frame Forces', ...
                   'NumberTitle', 'off', ...
                   'WindowState', 'maximized');
    h_ax = axes('Parent', h_fig);
    hold(h_ax, 'on'); axis equal; grid on; box off;
    xlabel(h_ax, 'X (Flight Direction)'); ylabel(h_ax, 'Y (Side)'); zlabel(h_ax, 'Z (Up)');
    title(h_ax, 'Body Orientation and Aerodynamic Forces');
    view(h_ax, 45, 30);

    ax_lim = params.c * 3;
    xlim(h_ax, [-ax_lim ax_lim]);
    ylim(h_ax, [-ax_lim ax_lim]);
    zlim(h_ax, [-ax_lim ax_lim]);

    % --- Cylinder geometry ---
    [Z_ext, Y_ext, X_ext] = cylinder(params.b, 50);
    X_ext = params.CoG(1) - X_ext * params.c;
    h_cyl = surf(h_ax, X_ext, Y_ext, Z_ext, ...
                 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'k', 'FaceAlpha', 0.6);
    Axis_X = [-1.5*params.c, 0.5*params.c] + params.CoG(1);
    h_axis = plot3(h_ax, Axis_X, [0 0], [0 0], 'r--', 'LineWidth', 1.5);
    h_cog = plot3(h_ax, 0, 0, 0, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 10);

    % --- Force arrows ---
    h_grav = quiver3(h_ax, 0, 0, 0, 0, 0, -params.m*params.g/force_scale, 'g', 'LineWidth', 2, 'MaxHeadSize', 0.5);
    h_lift = quiver3(h_ax, 0, 0, 0, 0, 0, 0, 'r', 'LineWidth', 2, 'MaxHeadSize', 0.5);
    h_drag = quiver3(h_ax, 0, 0, 0, 0, 0, 0, 'c', 'LineWidth', 2, 'MaxHeadSize', 0.5);

    % --- Wind arrows (blue wall) ---
    R_grid = 2 * params.c;
    N_arrows = 5;
    Y_positions = linspace(-R_grid, R_grid, N_arrows);
    Z_positions = linspace(-R_grid, R_grid, N_arrows);
    wind_start_x = 3 * params.c;
    wind_end_x = 2.5 * params.c;
    wind_vector_x = wind_end_x - wind_start_x;

    for i = 1:N_arrows
        for j = 1:N_arrows
            quiver3(h_ax, wind_start_x, Y_positions(i), Z_positions(j), ...
                    wind_vector_x, 0, 0, 0, 'b', 'LineWidth', 1.2, 'MaxHeadSize', 0.4);
        end
    end

    % --- Legend ---
    legend([h_grav, h_lift, h_drag, h_cyl], ...
           {'Gravity (Const.)', 'Lift (θ-frame)', 'Drag (θ-frame)', 'Cylinder'}, ...
           'Location', 'best');

    % --- Time slider and θ-checkbox ---
    t_max = length(SIM_DATA.alpha);
    uicontrol('Style', 'text', 'String', 'Time index:', ...
              'Position', [30 70 100 20], 'HorizontalAlignment', 'left');
    h_time = uicontrol('Style', 'slider', 'Min', 1, 'Max', t_max, 'Value', 1, ...
                       'SliderStep', [1/(t_max-1) 10/(t_max-1)], ...
                       'Position', [130 70 300 20], 'Callback', @updateVisualization);
    h_time_val = uicontrol('Style', 'edit', 'String', '1', ...
                           'Position', [450 70 60 25], 'Callback', @updateTimeEdit);
    h_no_sdslip_angle_checkbox = uicontrol('Style', 'checkbox', 'String', 'Include θ rotation', ...
                                 'Value', 1, ...
                                 'Position', [30 40 200 20], ...
                                 'Callback', @updateVisualization);

    % --- Text display for angles ---
    h_text = uicontrol('Style', 'text', ...
                       'String', '', ...
                       'Position', [30 100 500 25], ...
                       'FontWeight', 'bold', ...
                       'HorizontalAlignment', 'left');

    % --- Initialize ---
    updateVisualization();

    % === Update visualization ===
    function updateVisualization(~, ~)
        idx = round(get(h_time, 'Value'));
        set(h_time_val, 'String', num2str(idx));

        % Angles
        alpha = SIM_DATA.alpha(idx);
        beta = SIM_DATA.beta(idx);
        no_sdslip_angle = SIM_DATA.no_sdslip_angle(idx);
        alpha_equivalent = SIM_DATA.alpha_equivalent(idx);

        include_no_sdslip_angle = get(h_no_sdslip_angle_checkbox, 'Value');

        % --- Rotation matrices ---
        Ry = [cos(alpha), 0, -sin(alpha);
              0, 1, 0;
              sin(alpha), 0, cos(alpha)];
        Rz = [cos(beta), sin(beta), 0;
              -sin(beta),  cos(beta), 0;
              0, 0, 1];
        if include_no_sdslip_angle
            Rx = [1, 0, 0;
                  0, cos(no_sdslip_angle), -sin(no_sdslip_angle);
                  0, sin(no_sdslip_angle), cos(no_sdslip_angle)];
            R_total = Rz * Ry * Rx;
        else
            R_total = Rz * Ry;
        end

        % --- Rotate cylinder ---
        coords = [X_ext(:)'; Y_ext(:)'; Z_ext(:)'];
        coords_rot = R_total * coords;
        set(h_cyl, 'XData', reshape(coords_rot(1,:), size(X_ext)), ...
                   'YData', reshape(coords_rot(2,:), size(Y_ext)), ...
                   'ZData', reshape(coords_rot(3,:), size(Z_ext)));

        % --- Rotate axis ---
        axis_coords = [Axis_X; 0*Axis_X; 0*Axis_X];
        axis_rot = R_total * axis_coords;
        set(h_axis, 'XData', axis_rot(1,:), 'YData', axis_rot(2,:), 'ZData', axis_rot(3,:));

        % --- Aerodynamic forces in θ frame ---
        CL = params.CL(alpha_equivalent);
        CD = params.CD(alpha_equivalent);
        qbar = 0.5 * params.rho * SIM_DATA.V(idx)^2;

        F_L = qbar * params.S * [0; 0; CL];
        F_D = qbar * params.S * [-CD; 0; 0];

        F_L_rot = R_total * F_L;
        F_D_rot = R_total * F_D;

        % --- Update force arrows ---
        set(h_lift, 'UData', F_L_rot(1)/force_scale, ...
                    'VData', F_L_rot(2)/force_scale, ...
                    'WData', F_L_rot(3)/force_scale);
        set(h_drag, 'UData', F_D_rot(1)/force_scale, ...
                    'VData', F_D_rot(2)/force_scale, ...
                    'WData', F_D_rot(3)/force_scale);

        % --- Update text display ---
        set(h_text, 'String', sprintf('t=%d | α=%.2f° | β=%.2f° | θ=%.2f° | α_{total}=%.2f°', ...
               idx, rad2deg(alpha), rad2deg(beta), rad2deg(no_sdslip_angle), rad2deg(alpha_equivalent)));
    end

    % === Time edit box callback ===
    function updateTimeEdit(~, ~)
        val = str2double(get(h_time_val, 'String'));
        if isnan(val) || val < 1 || val > t_max
            set(h_time_val, 'String', num2str(round(get(h_time, 'Value'))));
            return;
        end
        set(h_time, 'Value', round(val));
        updateVisualization();
    end
end
