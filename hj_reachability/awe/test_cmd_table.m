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

%% Test the chi_a command and gamma_a commmand drive towards the figure 8
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

use_value_fcn = true;

if use_value_fcn
    load('tables.mat')
end

s       = mod(-pi/6/a0, 2*pi);
sigma   = 45;
h_tau   = 250/h0;
Va      = 31/v0;
chi_a   = 1.0407/a0;
gamma_a = 0.5697/a0;
Ft      = 1.600; % in kN
v_reelout = 3; % constant reeloutspeed of 3m/s

Lem.a       = 120;
Lem.b       = 200;
Lem.phi0    = 1;

if LongLatState
    extraArgs.direction = direction;
    extraArgs.Lem = Lem;
    [long, lat] = getLongLat(s* a0, sigma, h0 * h_tau, extraArgs);
    initialState = [long/a0, lat/a0, h_tau, Va, chi_a, gamma_a, Ft]';
else
    initialState = [s, sigma, h_tau, Va, chi_a, gamma_a, Ft]';
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
targetDistanceArgs.visualizeVec = false;
targetDistanceArgs.visualizePath = true;
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
axes1.XLim = [35 200];

% Create textbox
annotation(targetDistanceArgs.fig_handle,'textbox',...
    [0.447997776130467 0.205128205128205 0.0924025203854706 0.027972027972028],...
    'String','Starting Position',...
    'Interpreter','latex',...
    'FontSize',15,...
    'FitBoxToText','off',...
    'EdgeColor','none');

iterations = t_end/ dtSmall;
distance = zeros(1, iterations);
traj = nan(7, iterations);

targetDistanceArgs.visualize    = false;
targetDistanceArgs.AX = axes1;
[distance_tmp, sol,p_C_W, p_kite_W, extraOuts] = sys.getTargetdistance(sys.x, targetDistanceArgs);
plot3(p_kite_W{1}, p_kite_W{2}, p_kite_W{3}, 'bo', 'MarkerFaceColor','b', 'MarkerSize', 20)
sys.x(5) = extraOuts.chi_cmd_a/sys.a0;
sys.x(6) = extraOuts.gamma_cmd_a/sys.a0;

length_curve = 0;

for i = 1:iterations
    if mod(i-1,5) == 0
        targetDistanceArgs.visualize = true;
        if i==1
            targetDistanceArgs.visualizeVec = true;
        else
            targetDistanceArgs.visualizeVec = false;
        end
    else
        targetDistanceArgs.visualize = false;
        targetDistanceArgs.visualizeVec = false;
    end
    [distance_tmp, sol,p_C_W, p_kite_W, extraOuts] = sys.getTargetdistance(sys.x, targetDistanceArgs);
    distance(i) = distance_tmp;
    traj(:,i) = sys.x;
    
    if LongLatState
        long = sys.x(1);
        lat = sys.x(2);
        h_tau = sys.x(3);
    else
        [long, lat, t_tau, t_rot_tau, t_W, t_rot_W] = getLongLat(sys.x(1), sys.x(2), h0 * h_tau, targetDistanceArgs);
    end
    [pos_W_x,pos_W_y,pos_W_z] = sph2cart(long,lat,h0*h_tau);
    pos_before = [pos_W_x,pos_W_y,pos_W_z];
    if use_value_fcn
        u = get_path_u(sys.x', ...
            alpha_options, mu_options, alpha_max, grid_min, grid_max, dx, mu_max, I_table, alpha_min, mu_min);
        d = [v_reelout; zeros(3,1)];
        sys.updateState(u, dtSmall, sys.x, d);
    else
        sys.updateState([0; 0], dtSmall, sys.x, [0;0;0;0]);
    end

    if LongLatState
        long = sys.x(1);
        lat = sys.x(2);
        h_tau = sys.x(3);
    else
        [long, lat, t_tau, t_rot_tau, t_W, t_rot_W] = getLongLat(sys.x(1), sys.x(2), h0 * h_tau, targetDistanceArgs);
    end
    [pos_W_x,pos_W_y,pos_W_z] = sph2cart(long,lat,h0*h_tau);   
    pos_after = [pos_W_x,pos_W_y,pos_W_z];
    
    difference = pos_after-pos_before;
    length_curve = length_curve + norm(difference);
    %quiver3(pos_before(1), pos_before(2), pos_before(3), difference(1), difference(2), difference(3));
    
%     if LongLatState && i > 3 && any(max(abs(traj(1:2,1:i-1)-sys.x(1:2))) < 0.001, 'all')
%         disp('curve complete')
%         traj(:,i+1:end) = [];
%         break
%     elseif ~LongLatState && i > 3 && abs(traj(1,1)-mod(sys.x(1), 2*pi)) < 0.001
%         disp('curve complete')
%         traj(:,i+1:end) = [];
%         break
%     end
    
%     Only for no wind
%     if abs(norm(pos_before-pos_after) - sys.x(4) * dtSmall) > 0.01
%         warning('we should be travelling a distance of about va')
%     end
end
view(targetDistanceArgs.AX,[85.9885774495375 43.2006589259743]);
