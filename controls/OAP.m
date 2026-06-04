%%% Off Axis Propeller

function [T_control] = OAP(t, q, T_mag, alpha, beta, delta, R_ib, sim)

    %eul = quat2eul(q(:)', 'ZYX');
    %phi = eul(3);
    phi = atan2(-R_ib(2,3), R_ib(3,3));

    % Use the physical gyroscopic phase offset of 90 degrees.
    gyro = pi/2;
    alpha = beta; % using error instead of direction
    
    T_cmd = (1 - cos(phi + pi + gyro - alpha)) / 2 * T_mag;
    T_control = [0; T_cmd; 0];

end

