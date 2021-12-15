% Copyright (C) 2021  Nikolaus Vertovec
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
% :Revision: 14-December-2021
% :Author: Nikolaus Vertovec (nikolaus.vertovec@eng.ox.ac.uk)

function dOpt = optDstb(obj, ~, y, deriv, dMode, ~)
% dOpt = optDstb(obj, t, y, deriv, ~, ~)

if nargin < 5
  dMode = 'max';
end

if ~(strcmp(dMode, 'max') || strcmp(dMode, 'min'))
  error('dMode must be ''max'' or ''min''!')
end

%% Calculate d1
if ~iscell(deriv)
  deriv = num2cell(deriv);
end
d1 = -sign(deriv{4}) .* obj.max_tether_diff_dot;

%% Calculate d2-d4 
if isempty(obj.curve_direction)
    direction = 1;
else
    direction = obj.curve_direction;
end

s        = y{1};
sigma    = y{2};
h_tau    = y{3};
extraArgs.direction = direction;
[long, lat, ~, ~, t_W, t_rot_W] = getLongLat(s, sigma, obj.h0 * h_tau, extraArgs);
d_wind_max = obj.d_wind_max;

d2 = -sign(sin(lat) .* cos(long) .* (deriv{1} .* 2*pi/945 .* t_W{1}./(direction * norm_cellVec(t_W))  + ...
                                     deriv{2} .* t_rot_W{1}./(direction * norm_cellVec(t_rot_W))) + ...
                       sin(long) .* (deriv{1} .* 2*pi/945 .* t_W{2}./(direction * norm_cellVec(t_W))  + ...
                                     deriv{2} .* t_rot_W{2}./(direction * norm_cellVec(t_rot_W))) + ...
           cos(lat) .* cos(long) .* (deriv{1} .* 2*pi/945 .* t_W{3}./(direction * norm_cellVec(t_W))  + ...
                                     deriv{2} .* t_rot_W{3}./(direction * norm_cellVec(t_rot_W)) - ...
                                     deriv{3})) * d_wind_max;                       

d3 = sign(-sin(lat) .* sin(long) .* (deriv{1} .* 2*pi/945 .* t_W{1}./(direction * norm_cellVec(t_W))  + ...
                                     deriv{2} .* t_rot_W{1}./(direction * norm_cellVec(t_rot_W))) + ...
                       cos(long) .* (deriv{1} .* 2*pi/945 .* t_W{2}./(direction * norm_cellVec(t_W))  + ...
                                     deriv{2} .* t_rot_W{2}./(direction * norm_cellVec(t_rot_W))) - ...
           cos(lat) .* sin(long) .* (deriv{1} .* 2*pi/945 .* t_W{3}./(direction * norm_cellVec(t_W))  + ...
                                     deriv{2} .* t_rot_W{3}./(direction * norm_cellVec(t_rot_W)) - ...
                                     deriv{3})) * d_wind_max; 

d4 = sign(cos(lat) .* (deriv{1} .* 2*pi/945 .* t_W{1}./(direction * norm_cellVec(t_W))  + ...
                                     deriv{2} .* t_rot_W{1}./(direction * norm_cellVec(t_rot_W))) - ...
          sin(lat) .* (deriv{1} .* 2*pi/945 .* t_W{3}./(direction * norm_cellVec(t_W))  + ...
                                     deriv{2} .* t_rot_W{3}./(direction * norm_cellVec(t_rot_W)) - ...
                                     deriv{3})) * d_wind_max;
                                 
                                 
% deriv{1} .* 2*pi/945 .* t_W{1}./(direction * norm_cellVec(t_W))  + ...
%     deriv{2} .* t_rot_W{1}./(direction * norm_cellVec(t_rot_W)) + ...
%     deriv{3} .* cos(lat) .* cos(long)) * d_wind_max;  
% 
% d3 = -sign(deriv{1} .* 2*pi/945 .* t_W{2}./(direction * norm_cellVec(t_W))  + ...
%     deriv{2} .* t_rot_W{2}./(direction * norm_cellVec(t_rot_W)) + ...
%     cos(lat) .* sin(long).* deriv{3}) * d_wind_max;
% 
% d4 = -sign(deriv{1} .* 2*pi/945 .* t_W{3}./(direction * norm_cellVec(t_W))  + ...
%     deriv{2} .* t_rot_W{3}./(direction * norm_cellVec(t_rot_W)) + ...
%     sin(lat).* deriv{3}) * d_wind_max;

dOpt = {d1; d2; d3; d4};
end