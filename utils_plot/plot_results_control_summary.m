function plot_results_control_summary(t, x, sim, SIM_DATA)
%PLOT_RESULTS_CONTROL_SUMMARY Plot selected control/report figures only.
%   Includes:
%   - Body x velocity and spin rate
%   - Angular error and target roll
%   - 3D trajectory with waypoints
%   - OAP torque vectors
%   - Nutation angle

    % Define the save path
    saveDir = 'figures_rep';
    if ~exist(saveDir, 'dir')
        mkdir(saveDir);
    end

    u = x(:,1);
    p = x(:,4);
    pos = x(:,11:13);

    %% Figure 4: body x velocity, spin rate, and propeller rate
    figure('Name','v_{b,x}, \omega_{b,x}, and propeller rate');

    yyaxis left
    h1 = plot(t, u, 'r', 'LineWidth', 1.5);
    hold on;
    h2 = plot(t, p/(2*pi), 'Color', [0 0.5 0], 'LineWidth', 1.5, 'LineStyle', '-');
    ylabel('v_{b,x} [m/s] and \omega_{b,x} [Hz]', 'FontSize', 15);
    ax = gca;
    ax.YColor = 'k';

    yyaxis right
    if isfield(SIM_DATA, 'prop_RPM_cmd')
        h3 = plot(SIM_DATA.t, abs(SIM_DATA.prop_RPM_cmd), 'b', 'LineWidth', 1.5);
        ylabel('Propeller rate [RPM]', 'FontSize', 15); ax.YColor = 'b';
        legend([h1 h2 h3], ...
            'v_{b,x} [m/s]', ...
            '\omega_{b,x} [Hz]', ...
            '|\omega_{prop}| [RPM]');
    else
        h3 = [];
        ylabel('Propeller rate [RPM]');
        legend([h1 h2], 'v_{b,x} [m/s]', '\omega_{b,x} [Hz]');
    end

    xlabel('Time [s]', 'FontSize', 15);
    grid on;
    title('Body x-velocity, spin rate, and propeller rate', 'FontSize', 20);

    %% Figure 7: angular error and target roll
    if isfield(SIM_DATA, 'delta') && isfield(SIM_DATA, 'alpha_target')
        figure('Name','Angular Error and Target Roll');
        plot(t, rad2deg(SIM_DATA.delta), 'r', 'LineWidth', 1.5);
        hold on;
        plot(t, rad2deg(SIM_DATA.alpha_target), 'b--', 'LineWidth', 1.5);
        if isfield(SIM_DATA, 'beta_target')
            plot(t, rad2deg(SIM_DATA.beta_target), 'g-.', 'LineWidth', 1.5);
            legend('\delta (error)', '\alpha_{target}', '\beta_{target}');
        else
            legend('\delta (error)', '\alpha_{target}');
        end
        hold off;
        xlabel('Time [s]');
        ylabel('Angle [deg]');
        grid on;
        title('Angular error \delta and target roll \alpha');
    end

    %% Figure 11: trajectory
    figure('Name','Trajectory');
    plot3(pos(:,1), pos(:,2), pos(:,3), 'r', 'LineWidth', 1.5);
    hold on;
    if isfield(sim.options, 'control') && isfield(sim.options.control, 'waypoints')
        wps = sim.options.control.waypoints;
        if size(wps, 1) ~= 3
            wps = wps';
        end

        plot3(wps(1,:), wps(2,:), wps(3,:), 'ko', ...
            'MarkerFaceColor', 'y', ...
            'MarkerSize', 6, ...
            'DisplayName', 'Waypoints');

        for k = 1:size(wps, 2)
            text(wps(1,k) + 0.8, wps(2,k) + 0.8, wps(3,k) - 0.2, ...
                sprintf('WP%d', k), ...
                'Color', 'k', ...
                'FontSize', 9, ...
                'FontWeight', 'bold', ...
                'HorizontalAlignment', 'left', ...
                'VerticalAlignment', 'middle', ...
                'Interpreter', 'none');
        end

        if isfield(sim.options.control, 'loop') && sim.options.control.loop
            wps_line = [wps, wps(:,1)];
        else
            wps_line = wps;
        end
        plot3(wps_line(1,:), wps_line(2,:), wps_line(3,:), 'k--', ...
            'LineWidth', 0.8, ...
            'DisplayName', 'Waypoint path');
    end
    hold off;
    xlabel('X [m]', 'FontSize', 15);
    ylabel('Y [m]', 'FontSize', 15);
    zlabel('Z [m]', 'FontSize', 15);
    grid on;
    axis equal;
    title('Trajectory in inertial frame', 'FontSize', 20);
    set(gca, 'YDir', 'reverse', 'ZDir', 'reverse');

    %% Figure 20: OAP torque vectors
    if isfield(SIM_DATA, 'T_mag') && isfield(SIM_DATA, 'T_control')
        hTorque = figure('Name', 'OAP Torque Vectors'); % ASSIGN HANDLE HERE
        plot(SIM_DATA.t, SIM_DATA.T_mag, 'b', 'LineWidth', 1.5);
        hold on;
        plot(SIM_DATA.t, SIM_DATA.T_mag_undamped, 'k', 'LineWidth', 1.5);
        %plot(SIM_DATA.t, SIM_DATA.T_control, 'r', 'LineWidth', 1.5);
        hold off;

        xlabel('Time [s]', 'FontSize', 15);
        ylabel('Torque magnitude [Nm]', 'FontSize', 15);
        legend('|T|','|T_{undamped}|', '|T_{control}|');
        grid on;
        title('Control Torque', 'FontSize', 20);

        t_zoom = [3.0 4.0];
        idx = SIM_DATA.t >= t_zoom(1) & SIM_DATA.t <= t_zoom(2);
        if any(idx)
            axes('Position', [0.58 0.50 0.30 0.30]);
            box on;
            hold on;
            plot(SIM_DATA.t(idx), SIM_DATA.T_control(idx), 'r', 'LineWidth', 1.0);
            plot(SIM_DATA.t(idx), SIM_DATA.T_mag(idx), 'b', 'LineWidth', 1.0);
            plot(SIM_DATA.t(idx), SIM_DATA.T_mag_undamped(idx), 'k', 'LineWidth', 1.0);
            hold off;
            grid on;
            xlim(t_zoom);

            y_zoom = [
                SIM_DATA.T_mag(idx);
                SIM_DATA.T_mag_undamped(idx);
                SIM_DATA.T_control(idx)
            ];
            y_pad = 0.05 * max(range(y_zoom), 1e-6);
            ylim([min(y_zoom) - y_pad, max(y_zoom) + y_pad]);
            title('Zoom');
            set(gca, 'FontSize', 8);
            legend('|T_{control}|');
        end
        % Save the figure
        exportgraphics(hTorque, fullfile(saveDir, 'OAP_Torque_Vectors.png'), 'Resolution', 300);
    end

        %% Figure 20: OAP torque vectors
    if isfield(SIM_DATA, 'T_mag') && isfield(SIM_DATA, 'T_control')
        figure('Name', 'OAP Torque Vectors');
        plot(SIM_DATA.t, SIM_DATA.T_mag, 'b', 'LineWidth', 1.5);
        hold on;
        plot(SIM_DATA.t, SIM_DATA.T_mag_undamped, 'k', 'LineWidth', 1.5);
        %plot(SIM_DATA.t, SIM_DATA.T_control, 'r', 'LineWidth', 1.5);
        hold off;

        xlabel('Time [s]', 'FontSize', 15);
        ylabel('Torque magnitude [Nm]', 'FontSize', 15);
        legend('|T|','|T_{undamped}|', '|T_{control}|');
        grid on;
        title('Control Torque', 'FontSize', 20);
    end

    %% Figure 25: nutation angle
    plot_nutation_angle_only(t, x, sim);

end


function plot_nutation_angle_only(t, x, sim)

    % Define the save path
    saveDir = 'figures_rep';
    if ~exist(saveDir, 'dir')
        mkdir(saveDir);
    end
    
    p = x(:,4);
    q = x(:,5);
    r = x(:,6);

    Ix = sim.prop.Xzylo.I(1,1);
    Iy = sim.prop.Xzylo.I(2,2);
    Iz = sim.prop.Xzylo.I(3,3);

    if sim.options.propeller_on
        Ix = Ix + sim.prop.Motor.I(1,1) + sim.prop.Prop.I(1,1);
        Iy = Iy + sim.prop.Motor.I(2,2) + sim.prop.Prop.I(2,2);
        Iz = Iz + sim.prop.Motor.I(3,3) + sim.prop.Prop.I(3,3);
    end

    Hx = Ix * p;
    Hy = Iy * q;
    Hz = Iz * r;
    H_mag = sqrt(Hx.^2 + Hy.^2 + Hz.^2);
    nutation_deg = rad2deg(acos(max(-1, min(1, Hx ./ H_mag))));
    
    %H_b = [Hx; Hy; Hz];
    %v_b = [x(1); x(2); x(3)];
    %v_b = v_b / norm(v_b);
    %H_b = H_b / norm(H_b);

    %nutation_deg = rad2deg(acos(max(-1, min(1, dot(H_b, v_b)))));

    hNutation = figure('Name', 'Nutation Angle');
    plot(t, nutation_deg, 'b', 'LineWidth', 1.5);
    xlabel('Time [s]', 'FontSize', 15);
    ylabel('Nutation angle [deg]', 'FontSize', 15);
    title('Nutation angle: \angle(H, central axis)', 'FontSize', 20);
    grid on;

    t_zoom = [3.0 4.0];
    idx = t >= t_zoom(1) & t <= t_zoom(2);
    if any(idx)
        axes('Position', [0.58 0.55 0.30 0.30]);
        box on;
        plot(t(idx), nutation_deg(idx), 'b', 'LineWidth', 1.0);
        grid on;
        xlim(t_zoom);
        y_pad = 0.05 * max(range(nutation_deg(idx)), 1e-6);
        ylim([min(nutation_deg(idx)) - y_pad, max(nutation_deg(idx)) + y_pad]);
        title('Zoom');
        set(gca, 'FontSize', 8);
        
    end
    exportgraphics(hNutation, fullfile(saveDir, 'Nutation_Angle.png'), 'Resolution', 300);

    figure('Name', 'Nutation Angle');
    plot(t, nutation_deg, 'b', 'LineWidth', 1.5);
    xlabel('Time [s]', 'FontSize', 15);
    ylabel('Nutation angle [deg]', 'FontSize', 15);
    title('Nutation angle: \angle(H, central axis)', 'FontSize', 20);
    grid on;
    
end
