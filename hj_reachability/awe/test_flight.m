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

clear;
load('tables_HJB.mat')
%% Normalisation parameter
h0 = 1;
v0 = 1;
a0 = 1;
LongLatState = false;
direction = -1;

s       = -6*pi/6/a0;
sigma   = -15;
h_tau   = 450/h0;
Va      = 20/v0;
chi_a   = 1.0407/a0;
gamma_a = 0.5697/a0;
tether_diff = 2.5e-3;

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
sys.Ft_set = 1.652e+03;
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
targetDistanceArgs.visualizePath = false;
targetDistanceArgs.fig_num      = 1;
targetDistanceArgs.delta0       = 0.15;
targetDistanceArgs.direction    = direction;
 
sys.xhist = [];
dtSmall = 0.001;
t_end = 0.5;

targetDistanceArgs.fig_handle = figure(targetDistanceArgs.fig_num); 
% Create axes
axes1 = axes('Parent',targetDistanceArgs.fig_handle);
axis off;
hold(axes1,'on');
view(axes1,[85.9885774495375 43.2006589259743]);
xlabel(axes1, 'X','Visible','off')
ylabel(axes1, 'Y','Visible','off')
zlabel(axes1, 'Z','Visible','off')

figure2 = figure('Name', 'avoidance Maneuver'); 
axes2 = axes('Parent',figure2);
ylabel(axes2, 'Tether Force [kN]','Interpreter','latex')
xlabel(axes2, 'Time [s]','Interpreter','latex')
set(axes2,'FontSize',15,'GridAlpha',0.6,'GridLineStyle','--',...
        'MinorGridAlpha',0.2,'TickLabelInterpreter','latex');
hold(axes2,'on');

VS_1 = (grid_min(1) : dx(1) : grid_max(1))';
VS_2 = (grid_min(2) : dx(2) : grid_max(2))';
VS_3 = (grid_min(3) : dx(3) : grid_max(3))';
VS_4 = (grid_min(4) : dx(4) : grid_max(4))';
VS_5 = (grid_min(5) : dx(5) : grid_max(5))';
VS_6 = (grid_min(6) : dx(6) : grid_max(6))';
VS_7 = (grid_min(7) : dx(7) : grid_max(7))';

iterations = t_end/ dtSmall;
distance = zeros(1, iterations);
traj = nan(7, iterations);
Ft = nan(1, iterations);
%% Run
for counter = 1:20
    ft_norm = 0;
    while ft_norm < 1200 ||  ft_norm >= 1800
        initialState = grid_min + rand(7,1) .* (grid_max - grid_min);
        sys.x = initialState; 
        sys.x(1) = mod(sys.x(1), 2*pi);
        targetDistanceArgs.visualize    = false;
        targetDistanceArgs.AX = axes1;
        [distance_tmp, sol,p_C_W, p_kite_W, extraOuts] = sys.getTargetdistance(sys.x, targetDistanceArgs);
        %plot3(p_kite_W{1}, p_kite_W{2}, p_kite_W{3}, 'bo', 'MarkerFaceColor','b', 'MarkerSize', 20)
        sys.x(5) = extraOuts.chi_cmd_a/sys.a0;
        sys.x(6) = extraOuts.gamma_cmd_a/sys.a0;
        d1 = interpn(VS_1, VS_2, VS_3, VS_4, VS_5, VS_6, VS_7, d1_table, ...
                        sys.x(1), sys.x(2), sys.x(3), sys.x(4), sys.x(5), sys.x(6), sys.x(7), 'nearest');
        d2 = interpn(VS_1, VS_2, VS_3, VS_4, VS_5, VS_6, VS_7, d2_table, ...
                        sys.x(1), sys.x(2), sys.x(3), sys.x(4), sys.x(5), sys.x(6), sys.x(7), 'nearest');
        d3 = interpn(VS_1, VS_2, VS_3, VS_4, VS_5, VS_6, VS_7, d3_table, ...
                        sys.x(1), sys.x(2), sys.x(3), sys.x(4), sys.x(5), sys.x(6), sys.x(7), 'nearest');
        d4 = interpn(VS_1, VS_2, VS_3, VS_4, VS_5, VS_6, VS_7, d4_table, ...
                        sys.x(1), sys.x(2), sys.x(3), sys.x(4), sys.x(5), sys.x(6), sys.x(7), 'nearest');
        d = [d1;d2;d3;d4];
        [~, ~, ft] = sys.dynamics(0, sys.x, num2cell([0;0]), num2cell(d));
        ft_norm = norm_cellVec(ft);
    end
    
    length_curve = 0;
    
    for i = 1:iterations
        sys.x(1) = mod(sys.x(1), 2*pi);
        boundary_min = grid_min > sys.x;
        boundary_max = grid_max < sys.x;
        if any(boundary_min) || any(boundary_max)
            warning('out of bounds')
            break
        end
    
        if mod(i-1,5) == 0
            targetDistanceArgs.visualize = false;
            if i==iterations
                targetDistanceArgs.visualizePath = true;
                targetDistanceArgs.visualizeVec = false;
            else
                targetDistanceArgs.visualizeVec = false;
            end
        else
            targetDistanceArgs.visualize = false;
            targetDistanceArgs.visualizeVec = false;
        end
        if i==iterations
                targetDistanceArgs.visualizePath = true;
        end
        [distance_tmp, sol,p_C_W, p_kite_W, extraOuts] = sys.getTargetdistance(sys.x, targetDistanceArgs);
        distance(i) = distance_tmp;
    %     sys.x(5) = extraOuts.chi_cmd_a/sys.a0;
    %     sys.x(6) = extraOuts.gamma_cmd_a/sys.a0;
        traj(:,i) = sys.x;
        %ls_dot = max(min(dx(7), sys.max_reel_speed), -sys.max_reel_speed);
        
        if LongLatState
            long = sys.x(1);
            lat = sys.x(2);
            h_tau = sys.x(3);
        else
            [long, lat, t_tau, t_rot_tau, t_W, t_rot_W] = getLongLat(sys.x(1), sys.x(2), h0 * h_tau, targetDistanceArgs);
        end
        [pos_W_x,pos_W_y,pos_W_z] = sph2cart(long,lat,h0*h_tau);
        pos_before = [pos_W_x,pos_W_y,pos_W_z];
        
        I = interpn(VS_1, VS_2, VS_3, VS_4, VS_5, VS_6, VS_7, I_table, ...
                        sys.x(1), sys.x(2), sys.x(3), sys.x(4), sys.x(5), sys.x(6), sys.x(7), 'nearest');
        [I1,I2] = ind2sub([alpha_options, mu_options], I);
        alpha = alpha_min  + (I1-1)/(alpha_options-1) * (alpha_max-alpha_min);
        mu    = mu_min     + (I2-1)/(mu_options-1)    * (mu_max-mu_min);
        u = [alpha, mu];
        targetDistanceArgs.u = u;
    
        d1 = interpn(VS_1, VS_2, VS_3, VS_4, VS_5, VS_6, VS_7, d1_table, ...
                        sys.x(1), sys.x(2), sys.x(3), sys.x(4), sys.x(5), sys.x(6), sys.x(7), 'nearest');
        d2 = interpn(VS_1, VS_2, VS_3, VS_4, VS_5, VS_6, VS_7, d2_table, ...
                        sys.x(1), sys.x(2), sys.x(3), sys.x(4), sys.x(5), sys.x(6), sys.x(7), 'nearest');
        d3 = interpn(VS_1, VS_2, VS_3, VS_4, VS_5, VS_6, VS_7, d3_table, ...
                        sys.x(1), sys.x(2), sys.x(3), sys.x(4), sys.x(5), sys.x(6), sys.x(7), 'nearest');
        d4 = interpn(VS_1, VS_2, VS_3, VS_4, VS_5, VS_6, VS_7, d4_table, ...
                        sys.x(1), sys.x(2), sys.x(3), sys.x(4), sys.x(5), sys.x(6), sys.x(7), 'nearest');
        d = [d1;d2;d3;d4];
        [~, ~, ft] = sys.dynamics(0, sys.x, num2cell(u), num2cell(d));
        Ft(i) = norm_cellVec(ft);
        sys.updateState(u, dtSmall, sys.x, d);
       
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
    end
    plot(axes2, 0:dtSmall:t_end-dtSmall, smooth(Ft, 200),'LineWidth', 2)
end