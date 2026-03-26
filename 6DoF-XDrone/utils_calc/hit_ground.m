function [value_z, isterminal, direction] = hit_ground(~,x)
% Stops the simulation when z position <= 0
%
% Inputs:
%   t - current time
%   x - current state vector
%
% Outputs:
%   value - the function we monitor (stop when value = 0)
%   isterminal - 1 to stop the integrator
%   direction - 0 means all zero-crossings

value_z = x(13);   % z position in inertial frame; integration stops when value_z = 0
isterminal = 1;  % stop the integration
direction = 1;  % only detect decreasing zeros (falling through z=0)

end
