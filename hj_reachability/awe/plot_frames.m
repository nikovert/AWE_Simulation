% Copyright (C) 2022  Nikolaus Vertovec
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
% :Revision: 14-April-2022  
% :Author: Nikolaus Vertovec (nikolaus.vertovec@eng.ox.ac.uk)
%
%% Generate visualization of the standard reference frames
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
chi_a   = 0.5477/a0;
gamma_a = 0.1873/a0;
tether_diff = 3e-3;

alpha = 0.12;
mu_a = 0.3;

Lem.a       = 120;
Lem.b       = 200;
Lem.phi0    = 1;

extraArgs.Lem = Lem;
[long, lat] = getLongLat(s* a0, sigma, h0 * h_tau, extraArgs);

initialState = [s, sigma, h_tau, Va, chi_a, gamma_a, tether_diff]';

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

u = [alpha; mu_a];
beta = 0;
%% Define transfomation matrices
M_tauW = {-sin(lat).*cos(long), -sin(lat).*sin(long),  cos(lat);
          -sin(long)          , cos(long)           ,         0;
          -cos(lat).*cos(long), -cos(lat).*sin(long), -sin(lat)};
M_Wtau = transpose(M_tauW);

M_WO = {cos(sys.ENVMT.windDirection_rad), sin(sys.ENVMT.windDirection_rad), 0;
        sin(sys.ENVMT.windDirection_rad), -cos(sys.ENVMT.windDirection_rad), 0;
        0, 0, -1};
M_OW = transpose(M_WO);

M_AbarA = {1,0,0; 
           0, cos(mu_a), -sin(mu_a); 
           0, sin(mu_a), cos(mu_a)};
M_AbarO = {cos(chi_a) .* cos(gamma_a), sin(chi_a) .* cos(gamma_a), -sin(gamma_a); 
              -sin(chi_a),                 cos(chi_a),                             0; 
               cos(chi_a) .* sin(gamma_a), sin(chi_a) .* sin(gamma_a), cos(gamma_a)};
M_OAbar = transpose_cellMatrix(M_AbarO);

M_AB = {cos(alpha) .* cos(beta), sin(beta), sin(alpha) .* cos(beta); 
    -cos(alpha) .* sin(beta), cos(beta), -sin(alpha) .* sin(beta); 
    -sin(alpha), 0, cos(alpha)};

M_WAbar = mult_cellMatrix(M_WO, M_OAbar);

M_WA = mult_cellMatrix(M_WAbar, M_AbarA);

M_WB = mult_cellMatrix(M_WA, M_AB);
%% Set up figure
normalize    = false;
fig_num      = 1;

fig_handle = figure(fig_num); 
% Create axes
axes1 = axes('Parent',fig_handle);
axis off;
hold(axes1,'on');
view(axes1,[85.9885774495375 43.2006589259743]);
xlabel('X','Visible','off')
ylabel('Y','Visible','off')
zlabel('Z','Visible','off')

visualize([long(1), lat(1), h_tau(1), Va(1), chi_a(1), gamma_a(1)], u)
view(axes1,[89.100000018767 10.7999988816321]);
%%
[pos_G_W_x,pos_G_W_y,pos_G_W_z] = sph2cart(long,lat,h_tau);
p_kite_W = {pos_G_W_x;pos_G_W_y;pos_G_W_z};
scale = 30;
linewidth = 2;
arcfunc_args.scale = 0.92*scale;
%% Plot NED frame
NED_color = [0 0.4470 0.7410];
NED_ax_W_x = mult_cellMatrix(M_WO, [scale; 0; 0]);
NED_ax_W_y = mult_cellMatrix(M_WO, [0; scale; 0]);
NED_ax_W_z = mult_cellMatrix(M_WO, [0; 0; scale]);
quiver3(p_kite_W{1},p_kite_W{2},p_kite_W{3},NED_ax_W_x{1},NED_ax_W_x{2},NED_ax_W_x{3}, 'LineWidth', linewidth, 'Color', NED_color)
quiver3(p_kite_W{1},p_kite_W{2},p_kite_W{3},NED_ax_W_y{1},NED_ax_W_y{2},NED_ax_W_y{3}, 'LineWidth', linewidth, 'Color', NED_color)
quiver3(p_kite_W{1},p_kite_W{2},p_kite_W{3},NED_ax_W_z{1},NED_ax_W_z{2},NED_ax_W_z{3}, 'LineWidth', linewidth, 'Color', NED_color)

arcfunc_args.M = M_WO;
arcfunc_args.theta = 2*pi;
arc_points = arc_func(p_kite_W, 'xy', arcfunc_args);
plot3(arc_points(1,:), arc_points(2,:), arc_points(3,:), 'LineWidth', linewidth, 'Color', NED_color, 'LineStyle', '-.')

% Create textbox
annotation(fig_handle,'textbox',...
    [0.475999999999996 0.70238095238096 0.0329285714285711 0.054761904761905],...
    'String','$x_O$',...
    'Color', NED_color, ...
    'Interpreter','latex',...
    'FontSize',15,...
    'FitBoxToText','off',...
    'EdgeColor','none');

annotation(fig_handle,'textbox',...
    [0.808142857142852 0.626190476190486 0.0329285714285711 0.054761904761905],...
    'String','$y_O$',...
    'Color', NED_color, ...
    'Interpreter','latex',...
    'FontSize',15,...
    'FitBoxToText','off',...
    'EdgeColor','none');

annotation(fig_handle,'textbox',...
    [0.484928571428568 0.197619047619053 0.0329285714285711 0.0547619047619054],...
    'String','$z_O$',...
    'Color', NED_color, ...
    'Interpreter','latex',...
    'FontSize',15,...
    'FitBoxToText','off',...
    'EdgeColor','none');

%% Plot chi_a and gamma_a
arcfunc_args.M = M_WO;
arcfunc_args.theta = chi_a;
arc_points = arc_func(p_kite_W, 'xy', arcfunc_args);
plot3(arc_points(1,:), arc_points(2,:), arc_points(3,:), 'LineWidth', linewidth, 'Color', 'k', 'LineStyle', '-')
% Create textbox
annotation(fig_handle,'textbox',...
    [0.529571428571424 0.711904761904769 0.0329285714285711 0.054761904761905],...
    'String','$\chi_a$',...
    'Color', 'k', ...
    'Interpreter','latex',...
    'FontSize',20,...
    'FitBoxToText','off',...
    'EdgeColor','none');

arcfunc_args.M = M_WAbar;
arcfunc_args.theta = -gamma_a;
arc_points = arc_func(p_kite_W, 'zx', arcfunc_args);
plot3(arc_points(1,:), arc_points(2,:), arc_points(3,:), 'LineWidth', linewidth, 'Color', 'k', 'LineStyle', '-')
% Create textbox
annotation(fig_handle,'textbox',...
    [0.652785714285709 0.709523809523818 0.032928571428571 0.0547619047619051],...
    'String','$\gamma_a$',...
    'Color', 'k', ...
    'Interpreter','latex',...
    'FontSize',20,...
    'FitBoxToText','off',...
    'EdgeColor','none');
%% Plot A frame
A_color = [0.4660 0.6740 0.1880];
A_ax_W_x = mult_cellMatrix(M_WA, [scale; 0; 0]);
A_ax_W_y = mult_cellMatrix(M_WA, [0; scale; 0]);
A_ax_W_z = mult_cellMatrix(M_WA, [0; 0; scale]);
quiver3(p_kite_W{1},p_kite_W{2},p_kite_W{3},A_ax_W_x{1},A_ax_W_x{2},A_ax_W_x{3}, 'LineWidth', linewidth, 'Color', A_color)
quiver3(p_kite_W{1},p_kite_W{2},p_kite_W{3},A_ax_W_y{1},A_ax_W_y{2},A_ax_W_y{3}, 'LineWidth', linewidth, 'Color', A_color)
quiver3(p_kite_W{1},p_kite_W{2},p_kite_W{3},A_ax_W_z{1},A_ax_W_z{2},A_ax_W_z{3}, 'LineWidth', linewidth, 'Color', A_color)

arcfunc_args.M = M_WA;
arcfunc_args.theta = 2*pi;
arc_points = arc_func(p_kite_W, 'xy', arcfunc_args);
plot3(arc_points(1,:), arc_points(2,:), arc_points(3,:), 'LineWidth', linewidth, 'Color', A_color, 'LineStyle', '-.')

% Create textbox
annotation(fig_handle,'textbox',...
    [0.767071428571423 0.461904761904772 0.0329285714285711 0.0547619047619051],...
    'String','$y_A$',...
    'Color', A_color, ...
    'Interpreter','latex',...
    'FontSize',15,...
    'FitBoxToText','off',...
    'EdgeColor','none');

annotation(fig_handle,'textbox',...
    [0.425999999999997 0.242857142857148 0.0329285714285713 0.0547619047619054],...
    'String','$z_A$',...
    'Color', A_color, ...
    'Interpreter','latex',...
    'FontSize',15,...
    'FitBoxToText','off',...
    'EdgeColor','none');
%% Plot Abar frame
Abar_color = [0.8500 0.3250 0.0980];
Abar_ax_W_x = mult_cellMatrix(M_WAbar, [scale; 0; 0]);
Abar_ax_W_y = mult_cellMatrix(M_WAbar, [0; scale; 0]);
Abar_ax_W_z = mult_cellMatrix(M_WAbar, [0; 0; scale]);
quiver3(p_kite_W{1},p_kite_W{2},p_kite_W{3},Abar_ax_W_x{1},Abar_ax_W_x{2},Abar_ax_W_x{3}, 'LineWidth', linewidth, 'Color', Abar_color, 'LineStyle', ':')
quiver3(p_kite_W{1},p_kite_W{2},p_kite_W{3},Abar_ax_W_y{1},Abar_ax_W_y{2},Abar_ax_W_y{3}, 'LineWidth', linewidth, 'Color', Abar_color)
quiver3(p_kite_W{1},p_kite_W{2},p_kite_W{3},Abar_ax_W_z{1},Abar_ax_W_z{2},Abar_ax_W_z{3}, 'LineWidth', linewidth, 'Color', Abar_color)

arcfunc_args.M = M_WAbar;
arcfunc_args.theta = 2*pi;
arc_points = arc_func(p_kite_W, 'xy', arcfunc_args);
plot3(arc_points(1,:), arc_points(2,:), arc_points(3,:), 'LineWidth', linewidth, 'Color', Abar_color, 'LineStyle', '-.')

% Create textbox
annotation(fig_handle,'textbox',...
    [0.765285714285703 0.55714285714286 0.0311428571428658 0.0690476190476329],...
    'String','$y_{\overline{A}}$',...
    'Color', Abar_color, ...
    'Interpreter','latex',...
    'FontSize',15,...
    'FitBoxToText','off',...
    'EdgeColor','none');

annotation(fig_handle,'textbox',...
    [0.531357142857139 0.216666666666673 0.0329285714285711 0.0547619047619056],...
    'String','$z_{\overline{A}}$',...
    'Color', Abar_color, ...
    'Interpreter','latex',...
    'FontSize',15,...
    'FitBoxToText','off',...
    'EdgeColor','none');
%% Plot mu_a
arcfunc_args.M = M_WAbar;
arcfunc_args.theta = mu_a;
arc_points = arc_func(p_kite_W, 'zy', arcfunc_args);
plot3(arc_points(1,:), arc_points(2,:), arc_points(3,:), 'LineWidth', linewidth, 'Color', 'k', 'LineStyle', '-')
% Create textbox
annotation(fig_handle,'textbox',...
    [0.718857142857139 0.533333333333339 0.032928571428571 0.0547619047619053],...
    'String','$\mu_a$',...
    'Color', 'k', ...
    'Interpreter','latex',...
    'FontSize',20,...
    'FitBoxToText','off',...
    'EdgeColor','none');
%% Plot B frame
B_color = [0.6350 0.0780 0.1840];
B_ax_W_x = mult_cellMatrix(M_WB, [scale; 0; 0]);
B_ax_W_y = mult_cellMatrix(M_WB, [0; scale; 0]);
B_ax_W_z = mult_cellMatrix(M_WB, [0; 0; scale]);
quiver3(p_kite_W{1},p_kite_W{2},p_kite_W{3},B_ax_W_x{1},B_ax_W_x{2},B_ax_W_x{3}, 'LineWidth', linewidth, 'Color', B_color, 'LineStyle', '-')
quiver3(p_kite_W{1},p_kite_W{2},p_kite_W{3},B_ax_W_y{1},B_ax_W_y{2},B_ax_W_y{3}, 'LineWidth', linewidth, 'Color', B_color, 'LineStyle', ':')
quiver3(p_kite_W{1},p_kite_W{2},p_kite_W{3},B_ax_W_z{1},B_ax_W_z{2},B_ax_W_z{3}, 'LineWidth', linewidth, 'Color', B_color, 'LineStyle', '-')

arcfunc_args.M = M_WB;
arcfunc_args.theta = 2*pi;
arc_points = arc_func(p_kite_W, 'xy', arcfunc_args);
plot3(arc_points(1,:), arc_points(2,:), arc_points(3,:), 'LineWidth', linewidth, 'Color', B_color, 'LineStyle', '-.')
% Create textbox
annotation(fig_handle,'textbox',...
    [0.665285714285703 0.804761904761921 0.032928571428571 0.054761904761905],...
    'String','$x_B$',...
    'Color', B_color, ...
    'Interpreter','latex',...
    'FontSize',15,...
    'FitBoxToText','off',...
    'EdgeColor','none');

annotation(fig_handle,'textbox',...
    [0.470642857142852 0.259523809523817 0.0329285714285711 0.0547619047619059],...
    'String','$z_B$',...
    'Color', B_color, ...
    'Interpreter','latex',...
    'FontSize',15,...
    'FitBoxToText','off',...
    'EdgeColor','none');
%% Plot alpha
arcfunc_args.M = M_WB;
arcfunc_args.theta = -alpha;
arc_points = arc_func(p_kite_W, 'zx', arcfunc_args);
plot3(arc_points(1,:), arc_points(2,:), arc_points(3,:), 'LineWidth', linewidth, 'Color', 'k', 'LineStyle', '-')
% Create textbox
annotation(fig_handle,'textbox',...
    [0.654571428571424 0.773809523809532 0.032928571428571 0.0547619047619053],...
    'String','$\alpha_a$',...
    'Color', 'k', ...
    'Interpreter','latex',...
    'FontSize',20,...
    'FitBoxToText','off',...
    'EdgeColor','none');
%%
function [points] = arc_func(PoI, plane, extraArgs)
    if nargin < 3
        extraArgs = [];
    end
    
    if isfield(extraArgs, 'scale')
        scale = extraArgs.scale;
    else
        scale = 10;
    end

    if isfield(extraArgs, 'theta')
        theta = extraArgs.theta;
    else
        theta = 2*pi;
    end

    if isfield(extraArgs, 'nrpoints')
        nrpoints = extraArgs.nrpoints;
    else
        nrpoints = 500;
    end

    if isfield(extraArgs, 'M')
        M_frame = extraArgs.M;
        if iscell(M_frame)
            M_frame = cell2mat(M_frame);
        end
    else
        M_frame = eye(3);
    end
    if iscell(PoI)
        PoI = cell2mat(PoI);
    end
    switch plane
        case {'xy','yx'}
            M_rot = @(t) [cos(t), -sin(t), 0; ...
                     sin(t), cos(t), 0; ...
                     0, 0, 1];
            vec = [scale; 0; 0];
        case {'xz','zx'}
            M_rot = @(t) [cos(t), 0, sin(t); ...
                     0, 1, 0; ...
                     -sin(t), 0, cos(t)];
            vec = [scale; 0; 0];
        case {'yz','zy'}
            M_rot = @(t) [1, 0, 0; ...
                     0, cos(t), -sin(t); ...
                     0, sin(t), cos(t)];
            vec = [0; scale; 0];
        otherwise
            error('unknown plane')
    end
    points = zeros(3, nrpoints);
    for i = 1:nrpoints
        points(:,i) = PoI + M_frame * M_rot(i * theta/nrpoints) * vec;
    end
    
end

