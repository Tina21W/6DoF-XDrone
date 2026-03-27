function [T_control] = OAP(t, q)

    eul = quat2eul(q(:)', 'ZYX');
    phi = eul(3);

    T_control = ((sin(phi+pi) + 1) / 2)^2 * [0;0.01;0];

end
