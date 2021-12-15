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

%% Plot tether graphic
% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

%% Define and plot aircraft
pn = -185.9067; % inertial North position
pe = 173.7074; % inertial East position
pd = -200.0611; % inertial down position
   
vn = -8.8778;
ve = 14.8823;
vd = 3.1132;

phi = 0.6355;
theta = 0.2500;
psi = 0.6532;

scale = 10; 
[~, V, ~] = rndread('kitemill2.stl');
Vert = scale*V';
[p] = drawVehicleBody2(Vert,pn, pe, pd, phi, theta, psi,[],'normal');

xlabel('$$(m)$$', 'interpreter', 'latex')
ylabel('$$(m)$$', 'interpreter', 'latex')
zlabel('$$(m)$$', 'interpreter', 'latex')
set(gca, 'Fontsize' ,14);
grid on
box on;
axis equal; hold on;
set(gca,'TickLabelInterpreter','latex');

%% Define and plot tether
n_t_p = 5;
vel_scale = 3;

p_vec = [32.3117;28.2346;32.7669;64.3501;56.5633;65.7206;95.9552;85.1388;98.8793;126.9224;114.0706;132.3298;157.0283;143.5677;166.0718];
v_vec = [-1.2475;2.5519;0.3578;-2.6731;5.0506;0.9149;-3.8850;7.6860;1.1099;-5.5903;10.0585;1.9361;-7.0147;12.6197;2.2563];

v_vec = v_vec * vel_scale;
scatter3(0, 0, 0,'MarkerEdgeColor','k',...
        'MarkerFaceColor',[0 .75 .75])
    
for i = 0:n_t_p-1
    scatter3(p_vec(3*i+1), p_vec(3*i+2), p_vec(3*i+3),'MarkerEdgeColor','k',...
        'MarkerFaceColor',[0 .75 .75])
        quiver3(p_vec(3*i+1), p_vec(3*i+2), p_vec(3*i+3), v_vec(3*i+1), v_vec(3*i+2), v_vec(3*i+3), 'k')
end

%% Plot tether approximation
[long, lat, h_tau] = cart2sph(-pn, pe, -pd);
[long_dot, lat_dot, h_tau_dot] = vel_cart2sph(-pn, pe, -pd, vn, ve, vd);
p_vec_alt = zeros(size(p_vec));
v_vec_alt = zeros(size(p_vec));
for i = 1:3:n_t_p*3
   [pos_x,pos_y,pos_z] = sph2cart(long,lat,ceil(i/3) * h_tau/(n_t_p + 1));
   [vel_x,vel_y,vel_z] = vel_sph2cart(long,lat,ceil(i/3) * h_tau/(n_t_p + 1), long_dot,lat_dot,h_tau_dot);
   p_vec_alt(i)   = pos_x;
   p_vec_alt(i+1) = pos_y;
   p_vec_alt(i+2) = pos_z;
   
   v_vec_alt(i)   = vel_x;
   v_vec_alt(i+1) = vel_y;
   v_vec_alt(i+2) = vel_z;
    
   v_vec_alt(i:i+2) = v_vec_alt(i:i+2) * vel_scale;
   
   scatter3(p_vec_alt(i), p_vec_alt(i+1), p_vec_alt(i+2),'MarkerEdgeColor','k',...
        'MarkerFaceColor',[1 .0 .0])
    quiver3(p_vec_alt(i), p_vec_alt(i+1), p_vec_alt(i+2), v_vec_alt(i), v_vec_alt(i+1), v_vec_alt(i+2), 'r')
end

%% Plot lines between points
p_vec = [0;0;0;p_vec;-pn; pe; -pd];
p_vec_alt = [0;0;0;p_vec_alt;-pn; pe; -pd];
p_vec = reshape(p_vec, [3,n_t_p+2]);
p_vec_alt = reshape(p_vec_alt, [3,n_t_p+2]);

plot3(p_vec(1,:), p_vec(2,:), p_vec(3,:), 'k-');
plot3(p_vec_alt(1,:), p_vec_alt(2,:), p_vec_alt(3,:), ':r');

%% Extra Plot stuff
view(axes1,[-57.6170418167352 35.3859856017908]);
box(axes1,'on');
grid(axes1,'on');
axis(axes1,'tight');
hold(axes1,'off');
% Set the remaining axes properties
set(axes1,'DataAspectRatio',[1 1 1],'FontSize',14,'TickLabelInterpreter',...
    'latex');
% Create textbox
annotation(figure1,'textbox',...
    [0.567071428571426 0.111904761904765 0.0507857142857143 0.0452380952380955],...
    'String','$\mathbf{p}_{1}$',...
    'Interpreter','latex',...
    'FontSize',15,...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(figure1,'textbox',...
    [0.556357142857142 0.240476190476193 0.0507857142857143 0.0452380952380952],...
    'String','$\mathbf{p}_{2}$',...
    'Interpreter','latex',...
    'FontSize',15,...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(figure1,'textbox',...
    [0.542071428571427 0.400000000000003 0.0507857142857142 0.0452380952380955],...
    'String','$\mathbf{p}_{3}$',...
    'Interpreter','latex',...
    'FontSize',15,...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(figure1,'textbox',...
    [0.520642857142857 0.552380952380955 0.0507857142857142 0.0452380952380954],...
    'String','$\mathbf{p}_{4}$',...
    'Interpreter','latex',...
    'FontSize',15,...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(figure1,'textbox',...
    [0.504571428571428 0.704761904761907 0.0507857142857142 0.0452380952380951],...
    'String','$\mathbf{p}_{5}$',...
    'Interpreter','latex',...
    'FontSize',15,...
    'FitBoxToText','off',...
    'EdgeColor','none');
%% Functions
function pts = rotate(pts, phi, theta, psi)
    % From B 2 O
    pts = [ cos(psi)*cos(theta), cos(psi)*sin(phi)*sin(theta) - cos(phi)*sin(psi), sin(phi)*sin(psi) + cos(phi)*cos(psi)*sin(theta);
        cos(theta)*sin(psi), cos(phi)*cos(psi) + sin(phi)*sin(psi)*sin(theta), cos(phi)*sin(psi)*sin(theta) - cos(psi)*sin(phi);
        -sin(theta),                              cos(theta)*sin(phi),                              cos(phi)*cos(theta)] * pts;
end

function pts = translate(pts, pn, pe, pd)
    pts = pts + repmat([pn;pe;pd], 1, size(pts,2));
end

function [p] = drawVehicleBody2(V,pn, pe, pd, phi, theta, psi,p, mode)
    %V = rotate( V, phi+pi, theta, psi+pi );
    V = rotate( V, -phi-pi, -theta, psi+pi );
    % Stabilizer
    %V_stab = rotate(V_stab, phi, theta, psi); % body frame into NED frame
    V = translate(V, pn , pe, pd);
    M_EO = [-1, 0, 0; 0, 1, 0; 0, 0, -1];
    V = M_EO*V; % Transform from NED into E frame

    if isempty(p)
        [F, ~, C] = rndread('kitemill2.stl');
        col =  [ 57/255, 106/255, 177/255  ];
        p = patch('faces', F, 'vertices' ,V');
        set(p, 'facec', 'b');              % Set the face color (force it)
       % set(p, 'facec', 'flat');            % Set the face color flat
        set(p, 'FaceVertexCData', C);       % Set the color (from file)
        %set(p, 'facealpha',.4)             % Use for transparency
        set(p, 'EdgeColor','none');
        set(p, 'FaceLighting', 'none'); 
        set(p, 'FaceColor', col ); 
       % light  
        %('Position',[0 0 1]);% add a default light
        view(45,45);
        daspect([1 1 1])   
    else
        %set(p, 'Xdata', xstab, 'Ydata', ystab, 'Zdata', zstab);
        set(p, 'Vertices', V');
        %set(h_wing, 'Xdata', xwing, 'Ydata', ywing, 'Zdata', zwing);
        %set(h_plane, 'Xdata', xplane, 'Ydata', yplane, 'Zdata', zplane);
        drawnow
    end
end