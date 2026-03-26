clc, clear

% Input angles in degrees
beta_deg = 1;
alpha_deg = 1;
gamma_deg = 0;

% Convert to radians (Crucial Step!)
beta = deg2rad(beta_deg);
alpha = deg2rad(alpha_deg);
gamma = deg2rad(gamma_deg);

% Example quaternion
q = angle2quat(beta, alpha, gamma, 'ZYX');

% Normalize (important if q came from multiple operations)
q = q / norm(q);

% Extract scalar and vector parts
w = q(1);
v = q(2:4);

% Compute rotation angle
theta = 2 * acos(w);

% Compute rotation axis
if abs(sin(theta/2)) > 1e-6
    u = v / sin(theta/2);
else
    u = [0 0 0];  % arbitrary axis if rotation is nearly zero
end

mata = rad2deg(theta)
