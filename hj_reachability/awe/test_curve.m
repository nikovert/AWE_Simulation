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

%% Parameters
h_tau   = 250;
Lem.a       = 120/h_tau;
Lem.b       = 200/h_tau;
Lem.phi0    = 0;

long = @(s) Lem.b .* sin(s) ./( 1+(Lem.a/Lem.b .* cos(s)).^2 );
lat  = @(s) Lem.a .* sin(s).*cos(s) ./ ( 1+(Lem.a/Lem.b .* cos(s)).^2 );

%% Plot curve
s_range = 0:0.1:2*pi;
x = @(s) h_tau * cos(lat(s)) .* long(s);
y = @(s) h_tau * lat(s);

scatter(x(s_range),y(s_range),20)
hold on; grid on;

%% Derivatives
dlongds = @(s) ( Lem.b.^3 .* cos(s).*(2*Lem.a.^2-Lem.a.^2 .* cos(s).^2+Lem.b.^2)./(Lem.a.^2 .* cos(s).^2+Lem.b.^2).^2 );
dlatds  = @(s) ((cos(s).^2 .* (Lem.a.^3 .* Lem.b.^2+2*Lem.a .* Lem.b.^4) - Lem.a .* Lem.b.^4)./(Lem.a.^2 .* cos(s).^2+Lem.b.^2).^2 );

dxds = @(s) h_tau * (cos(lat(s)) .* dlongds(s) - sin(lat(s)) .* dlatds(s) .* long(s));
dyds = @(s) h_tau * dlatds(s);
quiver(x(s_range),y(s_range), dxds(s_range),dyds(s_range))

%% Test a point
extraArgs.Lem.a = 120;
extraArgs.Lem.b = 200;
extraArgs.Lem.phi0 = 0;

s_p = 2.5;

[long_p, lat_p, t_tau_p, t_rot_tau_p, t_W_p, t_rot_W_p] = getLongLat(s_p, 0, h_tau, extraArgs);
[x_p, y_p, z_p] = sph2cart(long_p, lat_p, h_tau);

disp('[0, x(s_p) y(s_p); x_p, y_p, z_p]:');
disp([0, x(s_p) y(s_p); x_p, y_p, z_p]);

disp('[0, dxds(s_p), dyds(s_p); cell2mat(t_W_p) * h_tau]:');
disp([0, dxds(s_p), dyds(s_p); cell2mat(t_W_p)' * h_tau]);

%% Test speed
v_k_tau = [25.5964 -11.7580 0]';

long_dot = (v_k_tau(2))./( abs(h_tau).*cos(lat(s_p)));
lat_dot = (v_k_tau(1))./abs(h_tau);
h_tau_dot = -(v_k_tau(3));

[vx, vy, vz] = vel_sph2cart(long_p, lat_p, h_tau, long_dot,lat_dot,h_tau_dot);

M_tauW = [-sin(lat_p) .* cos(long_p), -sin(lat_p) .* sin(long_p), cos(lat_p);
          -sin(long_p)            , cos(long_p)             , 0;
          -cos(lat_p) .* cos(long_p), -cos(lat_p) .* sin(long_p), -sin(lat_p) ];
      
v_k_W = M_tauW' * v_k_tau;
disp('[v_k_W; vx, vy, vz]:');
disp([v_k_W'; vx, vy, vz]);

%% Project to find s_dot
s_dot = ([dxds(s_p), dyds(s_p)] * [v_k_tau(1); v_k_tau(2)])/ (norm([dxds(s_p), dyds(s_p)]) * h_tau);

t_tau_p = cell2mat(t_tau_p);
s_dot_alt = (t_tau_p' * v_k_tau)/(norm(v_k_tau) * h_tau);


