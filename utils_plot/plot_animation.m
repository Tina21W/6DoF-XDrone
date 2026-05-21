function plot_animation(t, x, sim)
% PLOT_ANIMATION 6DoF Trajectory Animation (NED)
% Features:
%   - Precomputes rotation matrices for speed
%   - Quaternion normalization to prevent render warping
%   - GPU-Accelerated via hgtransform
%   - Simple 3D Cylinder with Heading Arrow

%% --- 1. Data Extraction and Setup ---
pos_i = x(:, 11:13); % Inertial Position (x, y, z)
quat  = x(:, 7:10);  % Quaternion
num_steps = length(t);
if num_steps == 0
    warning('Input state matrix (x) is empty. Cannot plot.');
    return;
end

%% --- 2. Interpolation to uniform time (50 FPS) ---
fps = 50;
t_uniform = linspace(t(1), t(end), round(fps*(t(end)-t(1))));

% Interpolate positions and quaternions
pos_i = interp1(t, pos_i, t_uniform, 'linear', 'extrap');
quat  = interp1(t, quat,  t_uniform, 'linear', 'extrap');

% Re-normalize quaternions to prevent hgtransform rendering warnings
quat_norms = sqrt(sum(quat.^2, 2)); 
quat = quat ./ quat_norms; 

t = t_uniform;
num_steps = length(t);

%% --- 3. Frame skipping option ---
frame_skip = 1;  

%% --- 4. Body visualization parameters ---
r_magnified = sim.options.radius_visualization;
chord = r_magnified * 1.5; % Length of the cylinder

% --- 3D Cylinder Grid (Base Frame at Origin) ---
num_circle_pts = 25; 
theta = linspace(0, 2*pi, num_circle_pts);
[Theta_mat, X_mat] = meshgrid(theta, [-chord/2, chord/2]);
Y_mat = r_magnified * cos(Theta_mat);
Z_mat = r_magnified * sin(Theta_mat);

% Blue/red lines along body z-axis (Base Frame at Origin)
blue_line_b_vector = [[0;0;0], [0;0;  r_magnified]];
red_line_b_vector  = [[0;0;0], [0;0; -r_magnified]];

% Central Arrow Parameters
arrow_start = -chord * 0.5;
arrow_end   = chord * 0.75;
arrow_head  = chord * 0.25; % Length of the arrow fins

%% --- 5. Precompute rotation matrices for speed ---
R_all = cell(num_steps,1);
for i = 1:num_steps
    R_all{i} = quatRotMat(quat(i,:));
end

%% --- 6. Figure Setup ---
fig_anim = figure('Name','6DoF Trajectory Animation','NumberTitle','off');
ax_anim = axes(fig_anim);

% Turn off hover toolbar safely (prevents hgtransform UI glitches)
try axtoolbar(ax_anim, 'off'); catch; end

hold(ax_anim,'on'); grid(ax_anim,'on');
xlabel(ax_anim,'North (m)'); ylabel(ax_anim,'East (m)'); zlabel(ax_anim,'Down (m)');
set(ax_anim,'YDir','reverse','ZDir','reverse'); % NED frame
view(ax_anim,45,20);
% Plot full trajectory
plot3(ax_anim,pos_i(:,1),pos_i(:,2),pos_i(:,3),'c--','DisplayName','Trajectory');

%% --- Waypoint visualization ---
if isfield(sim.options.control, 'waypoints')
    wps = sim.options.control.waypoints;  % 3 x N
    N_wp = size(wps, 2);
    colors_wp = lines(N_wp);  % distinct color per waypoint

    for k = 1:N_wp
        wx = wps(1,k); wy = wps(2,k); wz = wps(3,k);

        % Sphere marker
        plot3(ax_anim, wx, wy, wz, 'o', ...
            'MarkerSize', 12, ...
            'MarkerFaceColor', colors_wp(k,:), ...
            'MarkerEdgeColor', 'k', ...
            'LineWidth', 1.2, ...
            'DisplayName', sprintf('WP%d [%.0f,%.0f,%.0f]', k, wx, wy, wz));

        % Label offset slightly so it does not overlap the marker
        text(ax_anim, wx, wy, wz - sim.options.radius_visualization * 0.8, ...
            sprintf('  WP%d', k), ...
            'Color', colors_wp(k,:), ...
            'FontSize', 9, ...
            'FontWeight', 'bold', ...
            'HitTest', 'off');
    end

    % Connecting lines between waypoints (dashed, same color sequence)
    if isfield(sim.options.control, 'loop') && sim.options.control.loop
        idx_seq = [1:N_wp, 1];  % close the loop
    else
        idx_seq = 1:N_wp;
    end
    for k = 1:length(idx_seq)-1
        a = wps(:, idx_seq(k));
        b = wps(:, idx_seq(k+1));
        plot3(ax_anim, [a(1) b(1)], [a(2) b(2)], [a(3) b(3)], '--', ...
            'Color', [0.6 0.6 0.6], ...
            'LineWidth', 0.8, ...
            'HandleVisibility', 'off');
    end

    % Draw acceptance-radius circles (horizontal discs around each waypoint)
    if isfield(sim.options.control, 'waypoint_radius')
        r = sim.options.control.waypoint_radius;
        theta_circ = linspace(0, 2*pi, 40);
        for k = 1:N_wp
            cx = wps(1,k) + r*cos(theta_circ);
            cy = wps(2,k) + r*sin(theta_circ);
            cz = wps(3,k) * ones(1,40);
            plot3(ax_anim, cx, cy, cz, '-', ...
                'Color', [colors_wp(k,:), 0.35], ...
                'LineWidth', 0.8, ...
                'HandleVisibility', 'off');
        end
    end
end

%% --- 7. Compute extents (Uniform Box for Axis Equal) ---
x_min = min(pos_i(:,1)); x_max = max(pos_i(:,1));
y_min = min(pos_i(:,2)); y_max = max(pos_i(:,2));
z_min = min(pos_i(:,3)); z_max = max(pos_i(:,3));

x_span = (x_max - x_min) + 2*r_magnified;
y_span = (y_max - y_min) + 4*r_magnified;
z_span = (z_max - z_min) + 4*r_magnified;

x_center = (x_min + x_max)/2;
y_center = (y_min + y_max)/2;
z_center = (z_min + z_max)/2;

% max_span = max([x_span, y_span, z_span]);

% 1. Lock the 1:1:1 ratio first
axis(ax_anim, 'equal');

% 2. Force the limits second
set(ax_anim,'XLim', x_center + [-x_span/2 + arrow_start , x_span/2 + arrow_end], ...
            'YLim', y_center + [-y_span/2 + arrow_start, y_span/2 + arrow_end], ...
            'ZLim', z_center - chord/2 + [-z_span/2 , z_span/2]);

%% --- 8. Initial body graphics (HGTRANSFORM SETUP) ---
% Create a single transform group for the whole vehicle
tform_vehicle = hgtransform('Parent', ax_anim);

% Draw the Cylinder
surf(ax_anim, X_mat, Y_mat, Z_mat, ...
    'FaceColor',[0.6 0.6 0.8],'EdgeColor','k','FaceAlpha',0.8, ...
    'Parent', tform_vehicle, 'DisplayName','Body Cylinder');

% Draw the Z-Axis Indicator Lines
plot3(ax_anim, blue_line_b_vector(1,:), blue_line_b_vector(2,:), blue_line_b_vector(3,:), ...
    'b-','LineWidth',2, 'Parent', tform_vehicle, 'DisplayName','+Z Axis (Bottom)');

plot3(ax_anim, red_line_b_vector(1,:), red_line_b_vector(2,:), red_line_b_vector(3,:), ...
    'r-','LineWidth',2, 'Parent', tform_vehicle, 'DisplayName','-Z Axis (Top)');

% Draw the Central Arrow (Shaft)
plot3(ax_anim, [arrow_start, arrow_end], [0, 0], [0, 0], ...
    'g-', 'LineWidth', 2.5, 'Parent', tform_vehicle, 'DisplayName', 'Heading Axis (+X)');

% Draw the Central Arrow (Fins - XY Plane and XZ Plane)
plot3(ax_anim, [arrow_end - arrow_head, arrow_end, arrow_end - arrow_head], ...
               [arrow_head/2, 0, -arrow_head/2], [0, 0, 0], ...
    'g-', 'LineWidth', 2.5, 'Parent', tform_vehicle, 'HandleVisibility', 'off');

plot3(ax_anim, [arrow_end - arrow_head, arrow_end, arrow_end - arrow_head], ...
               [0, 0, 0], [arrow_head/2, 0, -arrow_head/2], ...
    'g-', 'LineWidth', 2.5, 'Parent', tform_vehicle, 'HandleVisibility', 'off');

title_h = title(ax_anim,sprintf('t = %.2f s',t(1)));
legend(ax_anim,'show','Location','bestoutside');

%% --- 9. Slider & Buttons ---
if num_steps > 1
    slider_h = uicontrol('Parent',fig_anim,'style','slider', ...
        'Min',1,'Max',num_steps,'Value',1,'Position',[250 20 400 20], ...
        'SliderStep',[1/(num_steps-1),10/(num_steps-1)]);
    addlistener(slider_h,'ContinuousValueChange',@(src,~) slider_callback(round(src.Value)));
else
    slider_h = uicontrol('Parent',fig_anim,'style','slider','Visible','off');
end
button_width = 80; button_height = 25; y_pos = 50;
uicontrol('Parent',fig_anim,'Style','pushbutton','String','▶ Play', ...
    'Position',[250 y_pos button_width button_height],'Callback',@(~,~) start_animation());
uicontrol('Parent',fig_anim,'Style','pushbutton','String','❚❚ Pause', ...
    'Position',[340 y_pos button_width button_height],'Callback',@(~,~) pause_animation());
uicontrol('Parent',fig_anim,'Style','pushbutton','String','⇤ Restart', ...
    'Position',[430 y_pos button_width button_height],'Callback',@(~,~) restart_animation());

%% --- 10. Timer ---
anim_timer = timer('ExecutionMode','fixedRate', ...
                   'Period',1/fps, ...
                   'TimerFcn',@(~,~) animate_step(), ...
                   'BusyMode','drop');
set(fig_anim,'DeleteFcn',@(~,~) delete_timer());

current_step = 1;
update_plot(current_step);

%% --- Nested Functions ---
    function update_plot(i)
        i = max(1, min(i, num_steps));
        P_i = pos_i(i,:)';
        R_i_b = R_all{i};
        
        % Build a 4x4 matrix combining Rotation and Translation
        M = eye(4);
        M(1:3, 1:3) = R_i_b; 
        M(1:3, 4)   = P_i;   
        
        % Move everything instantly via the GPU
        set(tform_vehicle, 'Matrix', M);
        
        set(title_h,'String',sprintf('t = %.2f s (step %d/%d)',t(i),i,num_steps));
        if num_steps>1, set(slider_h,'Value',i); end
        drawnow limitrate
    end

    function animate_step()
        if strcmp(anim_timer.Running,'on')
            if current_step >= num_steps
                pause_animation();
                set(title_h,'String','Animation complete.');
                return;
            end
            update_plot(current_step);
            current_step = current_step + frame_skip;
        end
    end
    function start_animation()
        if strcmp(anim_timer.Running,'off')
            if current_step >= num_steps, current_step = 1; end
            start(anim_timer);
        end
    end
    function pause_animation()
        if strcmp(anim_timer.Running,'on'), stop(anim_timer); end
    end
    function restart_animation()
        pause_animation();
        current_step = 1;
        update_plot(current_step);
    end
    function slider_callback(i)
        pause_animation();
        current_step = i;
        update_plot(current_step);
    end
    function delete_timer()
        if isvalid(anim_timer)
            stop(anim_timer);
            delete(anim_timer);
        end
    end
end