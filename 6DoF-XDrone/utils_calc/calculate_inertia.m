function results = calculate_inertia(m1, Ro1, Ri1, L1, x1, m2, Ro2, Ri2, L2, x2)
% CALCULATE_INERTIA Calculates the total Center of Mass (CoM) and the 
% 3x3 Inertia Tensor for a structure composed of two concentric annular wings (rings).
%
% This function assumes the structure is symmetric about the X-axis (y_i = z_i = 0),
% meaning the products of inertia (off-diagonal terms) will be zero.
%
% Inputs:
%   m1, m2 (double): Mass of Ring 1 and Ring 2 (kg).
%   Ro1, Ro2 (double): Outer Radius of Ring 1 and Ring 2 (m).
%   Ri1, Ri2 (double): Inner Radius of Ring 1 and Ring 2 (m).
%   L1, L2 (double): Length (chordwise width) of Ring 1 and Ring 2 (m).
%   x1, x2 (double): X-coordinate of the geometric center for Ring 1 and 2 (m).
%                    (Assumes CoM is at [x_i, 0, 0]).
%
% Output (struct):
%   results.M_total: Total mass (kg).
%   results.r_cg: Overall Center of Gravity position [x_cg; y_cg; z_cg] (m).
%   results.I_total: 3x3 Inertia Tensor about the overall CG (kg*m^2).

    % --- 1. Define Component Properties ---
    
    % Store masses and CoM positions
    masses = [m1; m2];
    r_i = [x1, 0.0, 0.0;
           x2, 0.0, 0.0]'; % Transposed to be 3x2 matrix [r1, r2]
    
    % Store dimensions
    Ro = [Ro1; Ro2];
    Ri = [Ri1; Ri2];
    L = [L1; L2];

    % --- 2. Calculate Total Mass and Overall Center of Gravity (CG) ---
    M_total = sum(masses);
    
    % Calculate the mass-weighted average position
    r_cg = (masses(1) * r_i(:, 1) + masses(2) * r_i(:, 2)) / M_total;

    % --- 3. Calculate Individual Inertia Tensors (I_i_geom) ---
    I_i_geoms = cell(2, 1);
    
    for i = 1:2
        m = masses(i);
        Ro_i = Ro(i);
        Ri_i = Ri(i);
        L_i = L(i);
        
        % Ixx: Roll Inertia (about the axis of symmetry)
        I_xx = 0.5 * m * (Ro_i^2 + Ri_i^2);
        
        % Iyy and Izz: Pitch/Yaw Inertia (about perpendicular axes)
        I_yy = 0.25 * m * (Ro_i^2 + Ri_i^2) + (1/12) * m * L_i^2;
        I_zz = I_yy;
        
        % Create the diagonal inertia tensor for component i about its own CoM
        I_i_geoms{i} = diag([I_xx, I_yy, I_zz]);
    end

    % --- 4. Apply Parallel Axis Theorem (PAT) to find Total Inertia (I_total) ---
    I_total = zeros(3, 3);
    E = eye(3); % Identity matrix

    for i = 1:2
        m = masses(i);
        I_i_geom = I_i_geoms{i};
        
        % Displacement vector from overall CG to component i's CoM (r_i')
        r_i_prime = r_i(:, i) - r_cg;
        
        % Calculate the PAT contribution matrix: m * ( (r' . r') * E - r' * r'^T )
        r_dot_r = dot(r_i_prime, r_i_prime); % (r' . r') (Dot product)
        r_rT = r_i_prime * r_i_prime';      % r' * r'^T (Matrix outer product)
        
        I_i_shift = m * (r_dot_r * E - r_rT);
        
        % Sum the shifted inertia to the total
        I_total = I_total + I_i_geom + I_i_shift;
    end

    % --- 5. Return Results ---
    results.M_total = M_total;
    results.r_cg = r_cg;
    results.I_total = I_total;
    
end