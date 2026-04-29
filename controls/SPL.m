%%% Swashplateless rotor

function [T_control] = SPL(t, q, T_matrix)

    % mode 1: Euler-roll version, breaks near the poles
    % mode 2: twist-quaternion version, breaks near 180 deg
    if nargin < 3
        mode = 2;
    end

    q = q(:);

    switch mode
        case 1
            eul = quat2eul(q(:)', 'ZYX');
            phi = eul(3);

            c_despun = cos(phi);
            s_despun = sin(phi);

            R_bds = [1 0 0
                     0 c_despun -s_despun
                     0 s_despun  c_despun];

        case 2
            q_twist = [q(1); q(2); 0; 0];

            if norm(q_twist) < 1e-12
                T_control = [0;0;0];
                return
            end

            q_twist = q_twist / norm(q_twist);
            R_bds = quatRotMat(q_twist);   % twist rotation about body x only

        otherwise
            error('SPL:InvalidMode', 'SPL mode must be 1 or 2.');
    end

    T_ds = T_matrix;           % T < 0.005 !!!
    T_control = R_bds' * T_ds;

end
