function plot_nutation_error(t, x, sim, p)
% PLOT_NUTATION_ERROR  Nutation angle of the X-zylo throughout the manoeuvre.
%
% Nutation angle is the angle between the angular momentum vector H and the
% body symmetry (spin) axis x_b = [1;0;0].  Computed entirely in the body
% frame, so no rotation matrices are needed:
%
%   H_b = I * omega_b  =  [Ix*p,  Iy*q,  Iz*r]
%   theta_n = acos( H_x / |H| )
%
% A second panel shows the instantaneous angular deviation of pitch and yaw
% from their slowly-varying means (Euler-angle amplitude) as a cross-check.

    % --- Extract body angular rates ---
    p = x(:,4);   % spin  (x-axis)
    q = x(:,5);   % pitch (y-axis)
    r = x(:,6);   % yaw   (z-axis)

    % --- Extract Euler angles from quaternion ---
    eul   = quat2eul(x(:,7:10), 'ZYX');
    psi   = eul(:,1);   % yaw
    theta = eul(:,2);   % pitch

    % --- Inertia components ---
    Ix = sim.prop.Xzylo.I(1,1);
    Iy = sim.prop.Xzylo.I(2,2);
    Iz = sim.prop.Xzylo.I(3,3);

    if sim.options.propeller_on
        Ix = Ix + sim.prop.Motor.I(1,1) + sim.prop.Prop.I(1,1);
        Iy = Iy + sim.prop.Motor.I(2,2) + sim.prop.Prop.I(2,2);
        Iz = Iz + sim.prop.Motor.I(3,3) + sim.prop.Prop.I(3,3);
    end

    % --- Angular momentum in body frame ---
    Hx = Ix * p;
    Hy = Iy * q;
    Hz = Iz * r;
    H_mag = sqrt(Hx.^2 + Hy.^2 + Hz.^2);

    % Nutation angle: angle between H and body x-axis (spin axis)
    nutation_deg = rad2deg(acos(Hx ./ H_mag));

    % --- Euler-angle oscillation amplitude (cross-check) ---
    % Remove slowly-varying mean with a moving average window (N/x of total time)
    N = length(t);
    win = max(round(N / 20), 1);
    theta_dev = theta - movmean(theta, win);
    psi_dev   = psi   - movmean(psi,   win);
    euler_amp_deg = rad2deg(sqrt(theta_dev.^2 + psi_dev.^2));

    % --- Gaussian smoothed euler angle spin ---
    spin_hz      = mean(abs(p)) / (2*pi);
    smooth_win   = max(0.5, 2 / max(spin_hz, 0.1));   % seconds
    theta_smooth = smoothdata(theta, 'gaussian', smooth_win, 'SamplePoints', t);
    psi_smooth   = smoothdata(psi,   'gaussian', smooth_win, 'SamplePoints', t);

    % % --- Angle between H and omega vectors ---
    % omega_mag = sqrt(p.^2 + q.^2 + r.^2);
    % H_dot_omega = Hx.*p + Hy.*q + Hz.*r;
    % H_omega_angle_deg = rad2deg(acos(H_dot_omega ./ (H_mag .* omega_mag)));

    % ========================= PLOT =========================================
    figure('Name', 'Nutation Error');

    % -- Panel 1: primary nutation angle --
    subplot(4,1,1);
    plot(t, nutation_deg, 'b', 'LineWidth', 1.5);
    xlabel('Time [s]');
    ylabel('Nutation angle [deg]');
    title('Nutation angle: \angle(H, spin axis)');
    grid on;

    % -- Panel 2: Euler angle oscillation amplitude --
    subplot(4,1,2);
    plot(t, euler_amp_deg, 'r', 'LineWidth', 1.5);
    xlabel('Time [s]');
    ylabel('Amplitude [deg]');
    title('Euler-angle oscillation amplitude  \surd(\theta_{dev}^2 + \psi_{dev}^2)  [cross-check]');
    grid on;
    
    % -- panel 3: Smoothed euler angles --
    subplot(4,1,3);
    plot(t, rad2deg(theta_smooth), 'b', t, rad2deg(psi_smooth), 'r', 'LineWidth', 2);
    xlabel('Time [s]'); ylabel('Euler angles [deg]');
    title('Euler angles (low-pass smoothed — nutation removed)');
    grid on;
    
    % -- panel 4: Roll rate --
    subplot(4,1,4)
    plot(t, p/(2*pi), 'r', 'LineWidth', 1.5);
    xlabel('Time [s]'); ylabel('Angular Velocities [Hz]');
    legend('p_{nr}'); 
    grid on; title('Angular Velocity p (spin rate)');
    
    % -- Separate figure for only Panel 1 --
    figure('Name', 'Nutation Angle');
    plot(t, nutation_deg, 'b', 'LineWidth', 1.5);    
    xlabel('Time [s]');
    ylabel('Nutation angle [deg]');
    title('Nutation angle: \angle(H, V_{b})');
    grid on;

    % --- Zoomed inset ---
    t_zoom = [2.0 3.0];
    idx = t >= t_zoom(1) & t <= t_zoom(2);

    axes('Position', [0.58 0.55 0.30 0.30]);
    box on;
    plot(t(idx), nutation_deg(idx), 'b', 'LineWidth', 1.2);
    grid on;
    xlim(t_zoom);
    ylim([min(nutation_deg(idx)), max(nutation_deg(idx))]);
    title('Zoom');
    set(gca, 'FontSize', 8);
end
