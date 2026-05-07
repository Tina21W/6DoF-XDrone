%%% Off Axis Propeller

function [T_control] = OAP(t, q, T_vector)

    eul = quat2eul(q(:)', 'ZYX');
    phi = eul(3);

    T_control = ((sin(phi+pi) + 1) / 2) * T_vector;  % experiment with different functions

end

