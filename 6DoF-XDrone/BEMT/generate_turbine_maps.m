function [V_GRID, RPM_GRID, Torque_Map, Drag_Map] = generate_turbine_maps(Geom)
% GENERATE_TURBINE_MAPS Creates simple aerodynamic look-up tables for a turbine.
% NO INDUCTION CALCULATIONS INCLUDED.
%
% Inputs:
%   Geom - Structure containing all geometry and design parameters
% Outputs:
%   V_GRID, RPM_GRID  - Meshgrids of the operating conditions
%   Torque_Map        - Aerodynamic torque lookup table [Nm]
%   Drag_Map          - Axial drag lookup table [N]
    %% ========================================================================
    %%                    PART 1: UNPACK GEOMETRY & SETUP
    %% ========================================================================
    R_hub = Geom.hub_radius;
    R_duct = Geom.span/2;
    c_base = Geom.hub_height - 2e-3; 
    B = Geom.no_blades;
    rho = Geom.rho;
    V_inf_design = Geom.design_speed;
    RPM_design = Geom.design_RPM;
    Omega_design = RPM_design * 2 * pi / 60;
     % Air Density [kg/m^3]
    % 1 = Constant Z (Chord Projection) -> Chord INCREASES with twist
    % 2 = Constant Chord (Twisted)      -> Chord is constant
    blade_mode = 1;
    rotation_position = 0;
    alpha_d = 0;
    N_span  = 50;           % Grid resolution
    x_span = linspace(R_hub, R_duct, N_span);
    
    % --- TWIST CALCULATION ---
    phi_design = atan((Omega_design * x_span) ./ V_inf_design);
    beta_phys = phi_design + alpha_d; 
    
    %% -------- IMPORT GEOMETRY (NACA0012) --------
    % 1. Get the folder where THIS function lives (BEMT)
    func_dir = fileparts(mfilename('fullpath')); 
    
    % 2. Go up one level to the main folder (6DoF-XDrone)
    main_dir = fileparts(func_dir);
    
    % 3. Go into the new data folder
    data_dir = fullfile(main_dir, 'aerodynamics_data');
    
    % Build the full path to the airfoil file
    airfoil_filename = fullfile(data_dir, 'NACA0012.dat');
    
    if isfile(airfoil_filename)
        raw_data = readmatrix(airfoil_filename, 'NumHeaderLines', 2);
        split_idx = find(diff(raw_data(:,1)) < -0.5, 1);
        if isempty(split_idx)
            error('Could not identify Upper/Lower split in airfoil file.');
        else
            Upper_X = raw_data(1:split_idx, 1); Upper_Y = raw_data(1:split_idx, 2);
            Lower_X = raw_data(split_idx+1:end, 1); Lower_Y = raw_data(split_idx+1:end, 2);
            X_segs = [Upper_X; 1.0; flipud(Lower_X)];
            Y_segs = [Upper_Y; 0.0; flipud(Lower_Y)];
            X_norm = X_segs'; Y_norm = Y_segs';
            idx_LE = 1; idx_TE = length(Upper_X) + 1; 
        end
    else
        error('Airfoil file %s not found!\nLooked in: %s', 'NACA0012.dat', airfoil_filename);
    end
    
    X_air_base = (X_norm - rotation_position) * c_base; 
    Y_air_base = Y_norm * c_base;
    
    %% -------- CALCULATE SCALED CHORD DISTRIBUTION --------
    c_dist = zeros(size(x_span)); 
    beta_root = beta_phys(1) + pi;
    Z_root = X_air_base .* cos(beta_root) + Y_air_base .* sin(beta_root);
    H_ref = abs(Z_root(idx_TE) - Z_root(idx_LE));
    
    for i = 1:N_span
        beta = beta_phys(i) + pi;
        Z_curr = X_air_base .* cos(beta) + Y_air_base .* sin(beta);
        H_curr = abs(Z_curr(idx_TE) - Z_curr(idx_LE));
        if blade_mode == 1 
            if H_curr > 1e-6; scale_factor = H_ref / H_curr; else; scale_factor = 1; end
        else
            scale_factor = 1;
        end
        c_dist(i) = c_base * scale_factor;
    end
    %% ========================================================================
    %%            PART 2: LOAD POLAR DATA
    %% ========================================================================
    % Build the full path to the CSV file using the data_dir
    csv_filename = fullfile(data_dir, 'xf-n0012-il-200000.csv');
    
    PolarData = readtable(csv_filename, 'NumHeaderLines', 10);
    [Alpha_data, uI] = unique(PolarData.Alpha);
    Cl_data = PolarData.Cl(uI); Cd_data = PolarData.Cd(uI);
    
    %% ========================================================================
    %%            PART 3: SIMPLE AERODYNAMIC MAPS (NO BEM)
    %% ========================================================================
    v_vec   = linspace(0, 50, 50);   
    rpm_vec = linspace(0, 5000, 50); 
    [V_GRID, RPM_GRID] = meshgrid(v_vec, rpm_vec);
    Torque_Map = zeros(size(V_GRID)); Drag_Map = zeros(size(V_GRID));
    
    for i = 1:size(V_GRID, 1)
        for j = 1:size(V_GRID, 2)
            V_inf = V_GRID(i,j); 
            Omega = RPM_GRID(i,j) * 2 * pi / 60;
            T_sum = 0; D_sum = 0;
            
            for k = 1:N_span-1
                r_local = (x_span(k) + x_span(k+1))/2; 
                dr = x_span(k+1) - x_span(k);
                c_local = (c_dist(k) + c_dist(k+1))/2;
                
                % SIMPLE INFLOW (No induction 'a' or 'ap')
                phi = atan((Omega * r_local) / V_inf);
                beta_local = (beta_phys(k) + beta_phys(k+1))/2;
                alpha_deg = rad2deg(beta_local - phi);
                
                Cl = interp1(Alpha_data, Cl_data, alpha_deg, 'linear', 'extrap');
                Cd = interp1(Alpha_data, Cd_data, alpha_deg, 'linear', 'extrap');
                
                % SIMPLE RELATIVE VELOCITY
                V_rel = sqrt(V_inf^2 + (Omega*r_local)^2);
                L = 0.5 * rho * V_rel^2 * c_local * Cl; 
                D = 0.5 * rho * V_rel^2 * c_local * Cd; 
                
                % TURBINE PHYSICS FORCES
                dF_axial_drag = L * sin(phi) + D * cos(phi);
                dF_driving_torque = L * cos(phi) - D * sin(phi);
                
                T_sum = T_sum + dF_driving_torque * r_local * dr;
                D_sum = D_sum + dF_axial_drag * dr;
            end
            Torque_Map(i,j) = T_sum * B; 
            Drag_Map(i,j) = D_sum * B;
        end
    end
    
    %% ========================================================================
    %%            PART 4: PLOTTING
    %% ========================================================================
    fprintf('\n--- Displaying Turbine Performance Maps ---\n');
    
    figure('Color','w','Name','Turbine Performance Maps','WindowStyle','docked'); 
    
    % Torque Map
    subplot(1,2,1); 
    [C,h] = contourf(V_GRID, RPM_GRID, Torque_Map, 40); 
    clabel(C,h); colorbar; colormap(jet); hold on;
    contour(V_GRID, RPM_GRID, Torque_Map, [0 0], 'w', 'LineWidth', 3);
    plot(V_inf_design, RPM_design, 'p', 'MarkerSize', 18, 'MarkerFaceColor',[1 0.84 0], 'MarkerEdgeColor','k');
    xlabel('V [m/s]'); ylabel('RPM'); title('TURBINE TORQUE [Nm]');
    
    % Drag Map
    subplot(1,2,2); 
    [C2,h2] = contourf(V_GRID, RPM_GRID, Drag_Map, 40); 
    clabel(C2,h2); colorbar; colormap(hot); hold on;
    plot(V_inf_design, RPM_design, 'p', 'MarkerSize', 18, 'MarkerFaceColor',[1 0.84 0], 'MarkerEdgeColor','k');
    xlabel('V [m/s]'); ylabel('RPM'); title('AXIAL DRAG [N]');
end