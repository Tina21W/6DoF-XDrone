%%% Off Axis Propeller

function [T_control] = OAP(t, q, T_matrix)

    eul = quat2eul(q(:)', 'ZYX');
    phi = eul(3);

    T_control = ((sin(phi+pi) + 1) / 2) * T_matrix;  % experiment with different functions

end

