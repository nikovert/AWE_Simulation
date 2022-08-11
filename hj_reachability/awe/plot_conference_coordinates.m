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

%% Add relevant files to path
addpath('../')
addToolbox;
conference_mode = true;
if conference_mode
    main_color =[1 1 1];
else
    main_color =[0 0 0];
end
%% Setup Parameters
h = 200;
s = 11*pi/17;
sigma = 30;
[long, lat, t_tau, t_rot_tau, t_W, t_rot_W] = getLongLat(s, sigma, h);
[x, y, z] = sph2cart(long, lat, h);

%% Figure handlers
% Create figure
figure1 = figure('WindowState','maximized');

% Create axes
axes1 = axes('Parent',figure1);
axis off
hold(axes1,'on');

%% Plot point
if ~conference_mode
    plot3(x,y,z,'ko','MarkerSize',10,'MarkerFaceColor',main_color)
end
axis([0, h, -h, h, 0, h])
xlabel('X','Visible','off')
ylabel('Y','Visible','off')
zlabel('Z','Visible','off')
hold on

%% Draw lines to point
syms r j
xr = r .* cos(lat) .* cos(long);
yr = r .* cos(lat) .* sin(long);
zr = r .* sin(lat);
fplot3(xr,yr,zr,[0 h], 'k', 'LineWidth', 2, 'Parent',axes1)
%fplot3(xr,yr,sym(0),[0 h],'k--', 'LineWidth', 2, 'Parent',axes1)
%fplot3(xr,yr,sym(z),[0 h],'k--', 'LineWidth', 2, 'Parent',axes1)
%fplot3(sym(x),sym(y),h*sin(j),[0 lat],'k--', 'LineWidth', 2, 'Parent',axes1)

%% Plot Sphere and Axis
[X,Y,Z] = sphere(55);
sphere_color = abs(main_color - [0.741176470588235 0.733333333333333 0.733333333333333]);
surf(X*h,Y*h,Z*h, 'FaceAlpha', 0.15, 'EdgeColor', sphere_color, ...
    'FaceColor','none','LineStyle', '--', 'Parent',axes1)
line([0 0], [0 0], [h 0],'Color','red','LineStyle','--', 'Parent',axes1,'LineWidth', 2)
line([0 0], [h -h], [0 0],'Color','red','LineStyle','--', 'Parent',axes1,'LineWidth', 2)
line([h 0], [0 0], [0 0],'Color','red','LineStyle','--', 'Parent',axes1,'LineWidth', 2)

%% Add long-lat description
syms j t
xa = h*cos(j)*cos(t);
ya = h*cos(j)*sin(t);
fsurf(xa,ya,0,[0 lat 0 long],'FaceColor','b','EdgeColor','none', 'FaceAlpha', 0.15, 'Parent',axes1)
syms u v
xp = u*cos(v)*cos(long);
yp = u*cos(v)*sin(long);
zp = u*sin(v);
fsurf(xp,yp,zp,[0 h 0 lat],'FaceColor','g','EdgeColor','none', 'FaceAlpha', 0.15, 'Parent',axes1)

%% Plot curve
% Parameters
h_tau       = h;
Lem.a       = 120/h_tau;
Lem.b       = 200/h_tau;
Lem.phi0    = 1;

M_WP = [cos(Lem.phi0),0, -sin(Lem.phi0);0, 1, 0; sin(Lem.phi0),0, cos(Lem.phi0)];

long_curve = @(s) Lem.b .* sin(s) ./( 1+(Lem.a/Lem.b .* cos(s)).^2 );
lat_curve  = @(s) Lem.a .* sin(s).*cos(s) ./ ( 1+(Lem.a/Lem.b .* cos(s)).^2 );

% Plot curve
num_samples = 200;
s_range = linspace(0,2*pi, num_samples);
[x_curve, y_curve, z_curve] = sph2cart(long_curve(s_range), lat_curve(s_range), h_tau);
curve_cart = M_WP * [x_curve; y_curve; z_curve];
curve_color = [0.501960784313725 0.501960784313725 0.501960784313725];
scatter3(curve_cart(1,:), curve_cart(2,:), curve_cart(3,:),20,'filled','MarkerEdgeColor','none',...
    'MarkerFaceColor',curve_color)

% Plot segment
if ~conference_mode
    num_samples = floor(1000 * s/2*pi);
    s_range = linspace(0,s, num_samples);
    [x_curve, y_curve, z_curve] = sph2cart(long_curve(s_range), lat_curve(s_range), h_tau);
    curve_cart = M_WP * [x_curve; y_curve; z_curve];
    segment_color = [0 0 1];
    scatter3(curve_cart(1,:), curve_cart(2,:), curve_cart(3,:),20,'filled','MarkerEdgeColor','none',...
        'MarkerFaceColor',segment_color)
end
%% Draw Tau-frame
M_tauW = [-sin(lat) .* cos(long), -sin(lat) .* sin(long), cos(lat);
          -sin(long)            , cos(long)             , 0;
          -cos(lat) .* cos(long), -cos(lat) .* sin(long), -sin(lat) ];
scale_factor = 40;
      
x_tau = scale_factor * M_tauW' * [1; 0; 0];
y_tau = scale_factor * M_tauW' * [0; 1; 0];
z_tau = scale_factor * M_tauW' * [0; 0; 1];
tau_color = [0.588235294117647 0.301960784313725 0.819607843137255];
quiver3(x,y,z,x_tau(1),x_tau(2),x_tau(3), 'Color',tau_color,'LineWidth',3)
quiver3(x,y,z,y_tau(1),y_tau(2),y_tau(3), 'Color',tau_color,'LineWidth',3)
quiver3(x,y,z,z_tau(1),z_tau(2),z_tau(3), 'Color',tau_color,'LineWidth',3)

%% Draw Gamma-frame
if ~conference_mode
    if iscell(t_W)
        t_W = scale_factor * cell2mat(t_W);
    else
        t_W = scale_factor * t_W;
    end
    if iscell(t_rot_W)
        t_rot_W = cell2mat(t_rot_W);
    end
    
    gamma_color = [0.929411764705882 0.690196078431373 0.129411764705882];
    quiver3(x,y,z,t_W(1),t_W(2),t_W(3), 'Color',gamma_color,'LineWidth',3)
    %quiver3(x,y,z,t_rot_W(1),t_rot_W(2),t_rot_W(3), 'Color','g','LineWidth',3)
    
    
    [long_noSigma, lat_noSigma, ~, ~, ~, t_rot_W_noSigma] = getLongLat(s, 0, h);
    if iscell(t_rot_W_noSigma)
        t_rot_W_noSigma = cell2mat(t_rot_W_noSigma);
    end
    
    [x_noSigma, y_noSigma, z_noSigma] = sph2cart(long_noSigma, lat_noSigma, h);
    [long_dot_noSigma, lat_dot_noSigma, ~] = vel_cart2sph(x_noSigma, y_noSigma, z_noSigma, ...
        t_rot_W_noSigma(1)/norm(t_rot_W_noSigma),t_rot_W_noSigma(2)/norm(t_rot_W_noSigma),t_rot_W_noSigma(3)/norm(t_rot_W_noSigma));
    X_initial_noSigma = [long_noSigma, long_dot_noSigma, lat_noSigma, lat_dot_noSigma];
    
    [long_dot, lat_dot, ~] = vel_cart2sph(x, y, z, ...
        t_rot_W(1)/norm(t_rot_W),t_rot_W(2)/norm(t_rot_W),t_rot_W(3)/norm(t_rot_W));
    X_initial = [long, long_dot, lat, lat_dot];
    
    % Plot fll curve
    options = odeset('MaxStep', 5);
    [t,X]=ode23s(@cir,[0,h*2*pi],X_initial_noSigma, options);
    u=X(:,1);
    v=X(:,3);
    [x_geo, y_geo, z_geo] = sph2cart(u, v, h);
    scatter3(x_geo, y_geo, z_geo, 20,gamma_color,'filled', 'MarkerEdgeColor','none','MarkerFaceAlpha',0.7)
    
    % Plot sigma segment
    sigma_color = [0.1 0.1 0.8];
    options = odeset('MaxStep', 1);
    [t,X]=ode23s(@cir,[0,sigma],X_initial_noSigma, options);
    u=X(:,1);
    v=X(:,3);
    [x_geo, y_geo, z_geo] = sph2cart(u, v, h);
    scatter3(x_geo, y_geo, z_geo, 20,sigma_color,'filled', 'MarkerEdgeColor','none','MarkerFaceAlpha',0.7)
    
    % Plot y-axis
    options = odeset('MaxStep', 1);
    [t,X]=ode23s(@cir,[0,40],X_initial, options);
    u=X(:,1);
    v=X(:,3);
    [x_geo, y_geo, z_geo] = sph2cart(u, v, h);
    scatter3(x_geo, y_geo, z_geo, 20,gamma_color,'filled', 'MarkerEdgeColor','none','MarkerFaceAlpha',0.7)
    
    % Plot tangent
    tangent_color = 'g';
    quiver3(x_noSigma,y_noSigma,z_noSigma,t_W(1),t_W(2),t_W(3), 'Color',tangent_color,'LineWidth',3, 'LineStyle','-.')
    quiver3(x_noSigma,y_noSigma,z_noSigma,-t_W(1),-t_W(2),-t_W(3), 'Color',tangent_color,'LineWidth',3, 'LineStyle','-.')
end
%% Add further annotation
view(axes1,[108.903466695156 14.0976794597227]);
hold(axes1,'off');
% Create textbox
if ~conference_mode
    annotation(figure1,'textbox',...
        [0.75061486482009 0.408102692564992 0.0279710289710291 0.0384615384615382],...
        'String','t',...
        'FontSize',30,...
        'FitBoxToText','off',...
        'Color', main_color, 'EdgeColor','none');
end

% Create textbox
if ~conference_mode
    annotation(figure1,'textbox',...
        [0.722142642597867 0.488158872340272 0.027971028971029 0.0384615384615382],...
        'String','\sigma',...
        'FontSize',30,...
        'FitBoxToText','off',...
        'Color', main_color, 'EdgeColor','none');
end
% Create textbox
annotation(figure1,'textbox',...
    [0.873260223193059 0.290316825160004 0.027971028971029 0.0384615384615383],...
    'String','y_W',...
    'FontSize',30,...
    'FitBoxToText','off',...
    'Color', main_color, 'EdgeColor','none');

% Create textbox
annotation(figure1,'textbox',...
    [0.419521062271062 0.161843597368698 0.027971028971029 0.0384615384615383],...
    'String','x_W',...
    'FontSize',30,...
    'FitBoxToText','off',...
    'Color', main_color, 'EdgeColor','none');

% Create textbox
annotation(figure1,'textbox',...
    [0.579781531486756 0.158102692564992 0.027971028971029 0.0384615384615383],...
    'String','\lambda',...
    'FontSize',30,...
    'FitBoxToText','off',...
    'Color', main_color, 'EdgeColor','none');

% Create textbox
annotation(figure1,'textbox',...
    [0.737597658146166 0.31931144326036 0.0279710289710289 0.0384615384615383],...
    'String','\phi',...
    'FontSize',30,...
    'FitBoxToText','off',...
    'Color', main_color, 'EdgeColor','none');

% Create textbox
annotation(figure1,'textbox',...
    [0.646948206022833 0.400714694336422 0.027971028971029 0.0384615384615383],...
    'String','h_\tau',...
    'FontSize',30,...
    'FitBoxToText','off',...
    'Color', main_color, 'EdgeColor','none');

% Create textbox
annotation(figure1,'textbox',...
    [0.617909572599873 0.882271978250031 0.027971028971029 0.0384615384615382],...
    'String','z_W',...
    'FontSize',30,...
    'FitBoxToText','off',...
    'Color', main_color, 'EdgeColor','none');

% Create textbox
if ~conference_mode
    annotation(figure1,'textbox',...
        [0.675614864820091 0.732540894812183 0.0279710289710291 0.0384615384615382],...
        'String','s',...
        'FontSize',30,...
        'FitBoxToText','off',...
        'Color', main_color, 'EdgeColor','none');
end

% Create textbox
annotation(figure1,'textbox',...
    [0.700185919801591 0.655637537998768 0.0279710289710291 0.0384615384615382],...
    'String','x_\tau',...
    'FontSize',30,...
    'FitBoxToText','off',...
    'Color', main_color, 'EdgeColor','none');

% Create textbox
annotation(figure1,'textbox',...
    [0.781728651530891 0.557820300594564 0.027971028971029 0.0384615384615383],...
    'String','y_\tau',...
    'FontSize',30,...
    'FitBoxToText','off',...
    'Color', main_color, 'EdgeColor','none');

% Create textbox
if ~conference_mode
    annotation(figure1,'textbox',...
        [0.732120113302949 0.651360997454411 0.0279710289710291 0.0384615384615382],...
        'String','y_\Gamma',...
        'FontSize',30,...
        'FitBoxToText','off',...
        'Color', main_color, 'EdgeColor','none');
    
    % Create textbox
    annotation(figure1,'textbox',...
        [0.670312300469019 0.550912580941655 0.0279710289710291 0.0384615384615383],...
        'String','x_\Gamma',...
        'FontSize',30,...
        'FitBoxToText','off',...
        'Color', main_color, 'EdgeColor','none');
    
    % Create textbox
    annotation(figure1,'textbox',...
        [0.712420420375645 0.423552130767238 0.027971028971029 0.0384615384615382],...
        'String','C',...
        'FontSize',30,...
        'FitBoxToText','off',...
        'Color', main_color, 'EdgeColor','none');
    
    % Create textbox
    annotation(figure1,'textbox',...
        [0.727003753708978 0.523271231890833 0.027971028971029 0.0384615384615382],...
        'String','K',...
        'FontSize',30,...
        'FitBoxToText','off',...
        'Color', main_color, 'EdgeColor','none');
end
if conference_mode
    set(gcf, 'color', '#002147');
end

%% Function definitions

function xp = cir(t,x)
    xp = zeros(4,1);
    xp(1) = x(2);
    xp(2) = 2*tan(x(3))*x(2)*x(4);
    xp(3) = x(4);
    xp(4) = -sin(x(3))*cos(x(3))*x(2)^2;
end