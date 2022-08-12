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
% :Revision: 12-August-2022
% :Author: Nikolaus Vertovec (nikolaus.vertovec@eng.ox.ac.uk)

function [A1_table, A3_table ,A2_table, b_table] = generate_CBF_LMI(g, data, tau, dynSys, extraArgs)
%GENERATE_LOOKUP_TABLE generates the looup table used during Simulation
if nargin < 5
    extraArgs = [];
end
   
if isfield(extraArgs, 'newgrid')  && size(data, g.dim+1) == 1
    data = regrid(data, g, extraArgs.newgrid);
    g = extraArgs.newgrid;
end

tauLength = length(tau);
clns = repmat({':'}, 1, g.dim);
derivFunc = @upwindFirstWENO5;

N = g.shape;
grid_min = g.min;
grid_max = g.max;    
dx       = g.dx;
for i = 1:g.dim
  if isfield(g, 'bdry') && isequal(g.bdry{i}, @addGhostPeriodic)
      N(i) = N(i) + 1;
      grid_max(i) = grid_max(i) + dx(i);
  end
  switch i
      case {1, 5, 6}
          dx(i) = dx(i) * dynSys.a0;
          grid_min(i) = grid_min(i) * dynSys.a0;
          grid_max(i) = grid_max(i) * dynSys.a0;
      case 3
          dx(i) = dx(i) * dynSys.h0;
          grid_min(i) = grid_min(i) * dynSys.h0;
          grid_max(i) = grid_max(i) * dynSys.h0;
      case 4
          dx(i) = dx(i) * dynSys.v0;
          grid_min(i) = grid_min(i) * dynSys.v0;
          grid_max(i) = grid_max(i) * dynSys.v0;
  end
          
end

A1_table    = zeros(N, 'single');
A2_table    = zeros(N, 'single');
A3_table    = zeros(N, 'single');
b_table     = zeros(N, 'single');

t = 1;
BRS_at_t = data(clns{:}, t);
Deriv = computeGradients(g, BRS_at_t, true(g.dim, 1), derivFunc);

q1 = Deriv{1}/dynSys.a0;
q2 = Deriv{2};
q3 = Deriv{3}/dynSys.h0;
q4 = Deriv{4}/dynSys.v0;
q5 = Deriv{5}/dynSys.a0;
q6 = Deriv{6}/dynSys.a0;
q7 = Deriv{7};
va      = g.xs{4} * dynSys.v0;
gamma_a = g.xs{6} * dynSys.a0;

A1 = q4 .* 1/dynSys.AIRCRAFT.mass * dynSys.h0/dynSys.v0^2;
A2 = q5 .* 1/dynSys.AIRCRAFT.mass .* 1./(va*dynSys.v0.*cos(gamma_a))* dynSys.h0/dynSys.v0;
A3 = q6 .* 1/dynSys.AIRCRAFT.mass .* -1./(va*dynSys.v0)* dynSys.h0/dynSys.v0;
[~, A1] = augmentPeriodicData(g, A1);
[~, A2] = augmentPeriodicData(g, A2);
[~, A3] = augmentPeriodicData(g, A3);
A1_table(clns{:})    = single(A1);
A2_table(clns{:})    = single(A2);
A3_table(clns{:})    = single(A3);

%% Compute b
d = dynSys.optDstb(t, g.xs, Deriv);
h_tau    = g.xs{3};
chi_a    = g.xs{5} * dynSys.a0;
tether_diff  = g.xs{7};

if isempty(dynSys.curve_direction)
    direction = 1;
else
    direction = dynSys.curve_direction;
end

if ~isempty(dynSys.Lem)
    extraArgs.Lem = dynSys.Lem;
end

extraArgs.direction = direction;

xi = dynSys.ENVMT.windDirection_rad;
if dynSys.v0 ~= 1
    warning('wind shear model does not use normalisation')
end

s        = g.xs{1}* dynSys.a0;
sigma    = g.xs{2};
[long, lat, t_tau, t_rot_tau, t_W, t_rot_W] = getLongLat(s, sigma, dynSys.h0 * h_tau, extraArgs);

v_w_O = wind_shear(h_tau .* sin(lat) ,xi, dynSys.base_windspeed);
M_tauW = {-sin(lat) .* cos(long), -sin(lat) .* sin(long), cos(lat);
          -sin(long)            , cos(long)             , 0;
          -cos(lat) .* cos(long), -cos(lat) .* sin(long), -sin(lat) };

M_AbarO = {cos(chi_a) .* cos(gamma_a), sin(chi_a) .* cos(gamma_a), -sin(gamma_a); 
          -sin(chi_a),                 cos(chi_a),                             0; 
           cos(chi_a) .* sin(gamma_a), sin(chi_a) .* sin(gamma_a), cos(gamma_a)};
M_OAbar = transpose_cellMatrix(M_AbarO);

M_OW = {cos(xi),  sin(xi),  0; 
        sin(xi), -cos(xi),  0; 
        0,              0, -1};
M_WO = M_OW; 

v_a_O = mult_cellMatrix(M_OAbar,{va;0;0});
v_k_O = add_cellMatrix(v_a_O, v_w_O); 
v_k_W = mult_cellMatrix(M_WO,v_k_O); 
v_k_W = add_cellMatrix(v_k_W, {d{2}; d{3}; d{4}});
v_k_tau = mult_cellMatrix(M_tauW,v_k_W);

s_dot     = cell2mat(mult_cellMatrix(transpose(t_W), mult_cellMatrix(v_k_W ,dynSys.v0))) ...
        ./ (direction * norm_cellVec(t_W));
s_dot = 2*pi/945 * s_dot;

sigma_dot = cell2mat(mult_cellMatrix(transpose(t_rot_W), mult_cellMatrix(v_k_W ,dynSys.v0))) ...
    ./ (norm_cellVec(t_rot_W));

long_dot = (v_k_tau{2})./( abs(h_tau).*cos(lat));
lat_dot = (v_k_tau{1})./abs(h_tau);
h_tau_dot = -(v_k_tau{3});

x_W{1} = long;
x_W{2} = lat;
x_W{3} = h_tau*dynSys.h0;
x_W{4} = long_dot * dynSys.v0/dynSys.h0;
x_W{5} = lat_dot * dynSys.v0/dynSys.h0;
x_W{6} = h_tau_dot*dynSys.v0;
x_W{7} = tether_diff;
[f_kite_G, ~] = dynSys.tether_forces(x_W);
F_t_W = mult_cellMatrix({-1}, f_kite_G);
f_kite_G_O = mult_cellMatrix(M_OW,F_t_W);
F_rest = mult_cellMatrix(M_AbarO, add_cellMatrix({0;0;dynSys.AIRCRAFT.mass * dynSys.ENVMT.g},f_kite_G_O));
va_dot = 1/dynSys.AIRCRAFT.mass * F_rest{1} * dynSys.h0/dynSys.v0^2;
chi_a_dot  = 1/dynSys.AIRCRAFT.mass * 1./(va*dynSys.v0.*cos(gamma_a)) .* F_rest{2} * dynSys.h0/dynSys.v0;
gamma_a_dot = 1/dynSys.AIRCRAFT.mass * -1./(va*dynSys.v0) .* F_rest{3} * dynSys.h0/dynSys.v0;
tether_diff_dot = d{1};

H_terms = q1 .* s_dot* dynSys.v0/dynSys.h0/dynSys.a0 + ...
    q2 .* sigma_dot* dynSys.v0/dynSys.h0 + ...
    q3 .* h_tau_dot + ...
    q4 .* va_dot + ...
    q5 .* chi_a_dot/dynSys.a0 + ...
    q6 .* gamma_a_dot/dynSys.a0 + ...
    q7 .* tether_diff_dot * dynSys.v0/dynSys.h0;

CBF_term = extraArgs.lambdaDiscount * BRS_at_t;

% Calculate termporal deriv
grid_full = createGrid([g.min; 0], [g.max; tau(end)], [N'; length(tau)], extraArgs.pdDims, true);
BRS_diff_t = computeGradients(grid_full, data, [false(g.dim, 1); 1], derivFunc);

b = H_terms + CBF_term + BRS_diff_t{end}(clns{:}, t);
[~, b] = augmentPeriodicData(g, b);
b_table(clns{:})    = single(b);

if isfield(extraArgs, 'saveFile') && extraArgs.saveFile
    save('tables.mat', 'A1_table', 'A2_table', 'A3_table', ...
        'b_table', 'grid_min', 'grid_max', 'dx', '-v7.3');
end
end




