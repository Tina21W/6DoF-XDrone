clc, clear

alpha = 10; %degrees
beta = -10; %degrees
combinedRotation_funct(alpha,beta)

function combinedRotation_funct(alpha_deg, beta_deg)
% plotCombinedRotation Calculates and plots the combined rotation of a sphere
%   This function takes pitch (alpha) and yaw (beta) in degrees,
%   calculates the equivalent single rotation axis (k) and angle (theta),
%   and visualizes the results on a sphere.

%% 1. Input Angles and Conversion
alpha = deg2rad(alpha_deg); % Pitch (Y-axis rotation)
beta = deg2rad(beta_deg);   % Yaw (Z-axis rotation)

fprintf('--- Input Angles ---\n');
fprintf('Pitch (alpha): %.2f degrees\n', alpha_deg);
fprintf('Yaw (beta):    %.2f degrees\n', beta_deg);
fprintf('--------------------\n');

%% 2. Rotation Matrices (Yaw followed by Pitch)
% R_z(beta) - Yaw around Z-axis
R_z = [cos(beta), -sin(beta), 0;
       sin(beta),  cos(beta), 0;
       0,          0,         1];

% R_y(alpha) - Pitch around Y-axis
R_y = [cos(alpha), 0, sin(alpha);
       0,          1, 0;
       -sin(alpha), 0, cos(alpha)];

% R_total = R_y * R_z (order matters!)
R_total = R_y * R_z;

%% 3. Calculate the Single Angle (theta)
% theta = arccos((Trace(R) - 1) / 2)
trace_R = trace(R_total);
theta = acos((trace_R - 1) / 2);
theta_deg = rad2deg(theta);

fprintf('--- Equivalent Single Rotation ---\n');
fprintf('Rotation Angle (theta): %.2f degrees\n', theta_deg);

%% 4. Calculate the Single Axis (k)
% k is proportional to [R32-R23; R13-R31; R21-R12]
kx_unnorm = R_total(3, 2) - R_total(2, 3);
ky_unnorm = R_total(1, 3) - R_total(3, 1);
kz_unnorm = R_total(2, 1) - R_total(1, 2);

% Normalize the vector k
k_unnorm = [kx_unnorm; ky_unnorm; kz_unnorm];
k = k_unnorm / norm(k_unnorm);

fprintf('Rotation Axis (k) Vector:\n');
fprintf('  kx: %.4f\n', k(1));
fprintf('  ky: %.4f\n', k(2));
fprintf('  kz: %.4f\n', k(3));
fprintf('----------------------------------\n');

%% 5. Visualization on a Sphere
figure('Name', 'Combined Rotation Visualization');
axis equal;
hold on;
grid on;
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
title(sprintf('Rotation: Pitch=%.0f°, Yaw=%.0f° (Equivalent Angle=%.1f°)', alpha_deg, beta_deg, theta_deg));

% Set plot limits
lim = 1.2;
xlim([-lim lim]);
ylim([-lim lim]);
zlim([-lim lim]);

% Draw Sphere (unit radius)
[sx, sy, sz] = sphere(30);
h = surf(sx, sy, sz, 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'FaceColor', [0.8 0.8 0.8]);

% Define initial unit vectors (reference frame of the sphere)
i_init = [1; 0; 0]; % Red (X)
j_init = [0; 1; 0]; % Green (Y)
k_init = [0; 0; 1]; % Blue (Z)

% Calculate final unit vectors
i_final = R_total * i_init;
j_final = R_total * j_init;
k_final = R_total * k_init;

% --- Plot Initial Frame (Dashed) ---
quiver3(0, 0, 0, i_init(1), i_init(2), i_init(3), 'r--', 'LineWidth', 1);
quiver3(0, 0, 0, j_init(1), j_init(2), j_init(3), 'g--', 'LineWidth', 1);
quiver3(0, 0, 0, k_init(1), k_init(2), k_init(3), 'b--', 'LineWidth', 1);
text(i_init(1)*1.05, i_init(2)*1.05, i_init(3)*1.05, 'X_0', 'Color', 'r', 'FontWeight', 'bold');

% --- Plot Final Frame (Solid) ---
quiver3(0, 0, 0, i_final(1), i_final(2), i_final(3), 'r', 'LineWidth', 2);
quiver3(0, 0, 0, j_final(1), j_final(2), j_final(3), 'g', 'LineWidth', 2);
quiver3(0, 0, 0, k_final(1), k_final(2), k_final(3), 'b', 'LineWidth', 2);
text(i_final(1)*1.05, i_final(2)*1.05, i_final(3)*1.05, 'X_F', 'Color', 'r', 'FontWeight', 'bold');

% --- Plot Rotation Axis (k) ---
quiver3(0, 0, 0, k(1), k(2), k(3), 'm', 'LineWidth', 3, 'MaxHeadSize', 0.5);
text(k(1)*1.05, k(2)*1.05, k(3)*1.05, 'Axis k', 'Color', 'm', 'FontWeight', 'bold');

view(3); % Set 3D view
view(30,25)
hold off;

end