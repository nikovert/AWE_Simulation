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

%% Generate Figure to explain Guidance Strategy
% Generates Figure OptimalCMD.png
%% Add relevant files to path
addpath('../')
addToolbox;

clear; close all

%% Normalisation parameter
h0 = 1;
v0 = 1;
a0 = 1;
LongLatState = false;
direction = -1;

s       = -pi/6/a0;
sigma   = 50;
h_tau   = 250/h0;
Va      = 31/v0;
chi_a   = -0.7477/a0;
gamma_a = 0.1873/a0;
tether_diff = 3e-3;

Lem.a       = 120;
Lem.b       = 200;
Lem.phi0    = 1;

if LongLatState
    extraArgs.direction = direction;
    extraArgs.Lem = Lem;
    [long, lat] = getLongLat(s* a0, sigma, h0 * h_tau, extraArgs);
    initialState = [long/a0, lat/a0, h_tau, Va, chi_a, gamma_a, tether_diff]';
else
    initialState = [s, sigma, h_tau, Va, chi_a, gamma_a, tether_diff]';
end
sys = AWE_3DOF(initialState);
sys.h0 = h0;    
sys.v0 = v0;
sys.a0 = a0;
sys.v_w_O = {-15.8031/v0   -0.0001/v0   -0.0000/v0}';
sys.Ft_set = 1.652199337711693e+03;
sys.v_ro_set = 0.01; % don't normalise
sys.F_T_max = 1.6649;
sys.skipTether = false;
sys.curve_direction = direction;
sys.ignoreTetherDiff = true;
sys.Lem     = Lem;
sys.LongLatState = LongLatState;
%% Define Target States
targetDistanceArgs.Lem = Lem;
targetDistanceArgs.distanceOnly = false;
targetDistanceArgs.headingOnly  = false;
targetDistanceArgs.normalize    = false;
%targetDistanceArgs.h_tau        = 250;
targetDistanceArgs.visualize = true;
targetDistanceArgs.visualizeVec = true;
targetDistanceArgs.visualizePath = false;
targetDistanceArgs.fig_num      = 1;
targetDistanceArgs.delta0       = 0.15;
targetDistanceArgs.direction    = direction;

sys.x = initialState;  
sys.xhist = [];
dtSmall = 0.01;
t_end = 6.5;

targetDistanceArgs.fig_handle = figure(targetDistanceArgs.fig_num); 
% Create axes
axes1 = axes('Parent',targetDistanceArgs.fig_handle);
axis off;
hold(axes1,'on');
view(axes1,[85.9885774495375 43.2006589259743]);
xlabel('X','Visible','off')
ylabel('Y','Visible','off')
zlabel('Z','Visible','off')

iterations = t_end/ dtSmall;
distance = zeros(1, iterations);
traj = nan(7, iterations);

% Plot part of the Figure 8
sol_spread = -1.:0.005:-0.25;
long_path = Lem.b./ h_tau .* sin(sol_spread) ./ ( 1+(Lem.a./Lem.b .* cos(sol_spread)).^2 );
lat_path  = Lem.a./ h_tau .* sin(sol_spread) .* cos(sol_spread) ./ ( 1+(Lem.a./Lem.b .* cos(sol_spread)).^2 ) ;
[X,Y,Z] = sph2cart(long_path,lat_path,h_tau);
M_WP = {cos(Lem.phi0),0, -sin(Lem.phi0);0, 1, 0; sin(Lem.phi0),0, cos(Lem.phi0)};
points = cell2mat(M_WP) * [X;Y;Z];
plot3(axes1, points(1,:),points(2,:),points(3,:), 'LineWidth', 3)

targetDistanceArgs.AX = axes1;
[distance_tmp, sol,p_C_W, p_kite_W, extraOuts] = sys.getTargetdistance(sys.x, targetDistanceArgs);
plot3(p_kite_W{1}, p_kite_W{2}, p_kite_W{3}, 'bo', 'MarkerFaceColor','b', 'MarkerSize', 2)
sys.x(5) = extraOuts.chi_cmd_a/sys.a0;
sys.x(6) = extraOuts.gamma_cmd_a/sys.a0;

view(targetDistanceArgs.AX,[85.9885774495375 43.2006589259743]);
