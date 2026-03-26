%%%%%%%%%%%%%%%%%%%% AERODYNAMICS DATA READING FILE %%%%%%%%%%%%%%%%%%%%%%%
% This files reads all the lift and drag coefficients
% It reads the CoP location as a function of AoA
% Afterwards, it smoothens the data and prepares it for the utils functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ============================================================
%  1. LIFT CURVE (C_L)
% =============================================================
Lift_graph = readmatrix('C_L_xzylo.txt');
alpha_deg_CL = Lift_graph(:,1);
C_L = Lift_graph(:,2);

alpha_rad_CL = deg2rad(alpha_deg_CL);
C_L_smooth = smoothdata(C_L, 'loess', 5);
C_L_interp = griddedInterpolant(alpha_rad_CL, C_L_smooth, 'pchip');

%% ============================================================
%  2. DRAG CURVE (C_D)
% =============================================================
Drag_graph = readmatrix('C_D_xzylo.txt');
alpha_deg_CD = Drag_graph(:,1);
C_D = Drag_graph(:,2);

alpha_rad_CD = deg2rad(alpha_deg_CD);
C_D_smooth = smoothdata(C_D, 'loess', 5);
C_D_interp = griddedInterpolant(alpha_rad_CD, C_D_smooth, 'pchip');

%% ============================================================
%  3. CENTER OF PRESSURE (COP)
% =============================================================

CoP_graph = readmatrix('COP_xzylo.txt');
alpha_deg_CoP = CoP_graph(:,1);
CoP = CoP_graph(:,2);

alpha_rad_CoP = deg2rad(alpha_deg_CoP);
CoP_smooth = smoothdata(CoP, 'loess', 6);
CoP_interp = griddedInterpolant(alpha_rad_CoP, CoP_smooth, 'pchip');

%% ============================================================
%  4. Pitching Moment about Aerodynamic Center
% =============================================================

% This part now is directly included in the body files, not needed now here
% anymore (like  the equations from above as well).