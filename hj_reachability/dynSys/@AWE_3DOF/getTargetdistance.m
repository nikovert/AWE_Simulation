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

function [distance, sol,p_C_W, p_kite_W, extraOuts] = getTargetdistance(obj, state, extraArgs)
%GETVTONLEMNISCATE Summary of this function goes here
%   Detailed explanation goes here

% state = long, lat, h_tau, Va, chi_a, gamma_a
if ~iscell(state)
    state = num2cell(state, 2)';
end

if nargin < 3
    extraArgs = [];
end

if ~isfield(extraArgs, 's_old')
    s_old = zeros(size(state{1}));
else
    if size(extraArgs.s_old) == size(state{1})
        s_old = extraArgs.s_old;
    else
        s_old = zeros(size(state{1}));
    end
end

if isempty(obj.curve_direction) || ~isfield(extraArgs, 'direction')
    direction = 1;
elseif ~isempty(obj.curve_direction)
    direction = obj.curve_direction;
    if isfield(extraArgs, 'direction') && direction ~= extraArgs.direction
        warning('direction does not match')
    else
        extraArgs.direction = direction;
    end
elseif isfield(extraArgs, 'direction')
    direction = extraArgs.direction;
end

if ~isfield(extraArgs, 'Lem') && isempty(obj.Lem)
    Lem.a = 120;
    Lem.b = 200;
    Lem.phi0 = 1;
     
%     sol_spread = 0:0.1:2*pi;
%     long_spread = Lem.b .* sin(sol_spread) ./ ( 1+(Lem.a./Lem.b .* cos(sol_spread)).^2 );
%     [~, long_index] = min(abs(long_spread./h_tau(:) - long(:)),[],2);
%     
%     s = sol_spread(long_index);
%     Lem.phi0 = lat - reshape(Lem.a .* sin(s) .* cos(s) ./ ( 1+(Lem.a./Lem.b .* cos(s)).^2 ), size(lat))./(0.5*h_tau);
else
    if ~isempty(obj.Lem)
        extraArgs.Lem = obj.Lem;
    end
    Lem = extraArgs.Lem;
end

if (isfield(extraArgs, 'LongLatState') && extraArgs.LongLatState) || obj.LongLatState
    long        = state{1} * obj.a0;
    lat         = state{2} * obj.a0;
else
    [long, lat] = getLongLat(state{1}* obj.a0, state{2}, obj.h0 * state{3}, extraArgs);
    if ~isfield(extraArgs, 's_old')
        s_old = state{1}* obj.a0;
    end
end

if ~isfield(extraArgs, 'h_tau')
    h_tau = state{3} * obj.h0;
else
    h_tau = extraArgs.h_tau;
end
Va          = state{4} * obj.v0;
chi_a       = state{5} * obj.a0;
gamma_a     = state{6} * obj.a0;

if abs(obj.ENVMT.windDirection_rad- pi) > 0.0001
    warning('These conditions have not been tested')
end

if ~isfield(extraArgs, 'normalize')
    normalize = true;
else
    normalize = extraArgs.normalize;
end

if ~isfield(extraArgs, 'checks')
    checks = 1;
    if (isfield(extraArgs, 'LongLatState') && extraArgs.LongLatState) || obj.LongLatState
        checks = 4;
    end
else
    checks = extraArgs.checks;
end

visualize_distance = isfield(extraArgs, 'visualize') && extraArgs.visualize && length(h_tau) == 1;


if isfield(extraArgs, 'fig_handle')
    set(0, 'CurrentFigure', extraArgs.fig_handle)
    view(60, 30)
elseif isfield(extraArgs, 'fig_num') 
    figure(extraArgs.fig_num)
    hold on
    view(60, 30)
end

if visualize_distance
    if isfield(extraArgs, 'u')
        u = extraArgs.u;
    else
        u = [0;0];
    end
    
    visualize([long(1), lat(1), h_tau(1), Va(1), chi_a(1), gamma_a(1)], u)
end
    
if ~isfield(extraArgs, 'delta0')
    delta0 = 0.15;
else
    delta0 = extraArgs.delta0;
end

% dx = obj.dynamics(0, state, {0,0});
% long_dot    = dx{1};
% lat_dot     = dx{2};
% h_tau_dot   = dx{3};
% [vel_k_kite_W_x,vel_k_kite_W_y,vel_k_kite_W_z] = vel_sph2cart(long,lat,h_tau, long_dot,lat_dot,h_tau_dot);
% vel_k_kite_W = {vel_k_kite_W_x;vel_k_kite_W_y;vel_k_kite_W_z};

[pos_G_W_x,pos_G_W_y,pos_G_W_z] = sph2cart(long,lat,1);
Lem.a = Lem.a ./ h_tau;
Lem.b = Lem.b ./ h_tau;

p_kite_W = {pos_G_W_x;pos_G_W_y;pos_G_W_z};
p_kite_W_norm = norm_cellVec(p_kite_W);
M_WP = {cos(Lem.phi0),0, -sin(Lem.phi0);0, 1, 0; sin(Lem.phi0),0, cos(Lem.phi0)};
pos_P = mult_cellMatrix(transpose(M_WP), p_kite_W); 

%% getTargetOnBoothNewton2 (p_kite_W is already normalised)
sol = cell(checks, 1);
p_C_P = cell(checks, 1);
nonconverged_points  = cell(checks, 1);
p_C_W  = cell(checks, 1);
p_C_W_norm = cell(checks, 1);
tmp = cell(checks, 1);
delta = cell(checks, 1);
for i = 1:checks
    [ sol{i},p_C_P{i}, nonconverged_points{i}] = getTargetOnBoothNewton2(Lem, pos_P, mod(s_old+(i-1)/checks*2*pi, 2*pi), h_tau, direction);
    p_C_W{i} = mult_cellMatrix(M_WP , p_C_P{i}); %LissajousFigure8Following Line 61

    p_C_W_norm{i} = norm_cellVec(p_C_W{i});
    % using the definition of the arc length on the unit sphere
    tmp{i} = (p_kite_W{1} .* p_C_W{i}{1} + p_kite_W{2} .* p_C_W{i}{2} + p_kite_W{3} .* p_C_W{i}{3})./(p_kite_W_norm .* p_C_W_norm{i});
    delta{i} = acos( min( max( tmp{i}, -1),1 ) );
end
p_kite_W = {p_kite_W{1} .* h_tau; p_kite_W{2} .* h_tau; p_kite_W{3} .* h_tau};

%%
t_P = cell(checks, 1);
t_W = cell(checks, 1);
for i = 1:checks
    t_P{i} = getBoothInfos2(sol{i},Lem, direction);
    t_W{i} = mult_cellMatrix(M_WP,t_P{i}); 
end

if visualize_distance && isfield(extraArgs, 'visualizePath') && extraArgs.visualizePath
    sol_spread = 0:0.01:2*pi;
    long_path = Lem.b .* sin(sol_spread) ./ ( 1+(Lem.a./Lem.b .* cos(sol_spread)).^2 );
    lat_path  = Lem.a .* sin(sol_spread) .* cos(sol_spread) ./ ( 1+(Lem.a./Lem.b .* cos(sol_spread)).^2 ) ;
    [X,Y,Z] = sph2cart(long_path,lat_path,h_tau);
    points = cell2mat(M_WP) * [X;Y;Z];

%     % Vectors along the points
%     t_range = nan(3,length(0:0.1:2*pi));
%     for i = 1:length(sol_spread)
%         t_P_tmp = getBoothInfos2(sol_spread(i),Lem, direction);
%         t_tmp = mult_cellMatrix(M_WP,t_P_tmp);
%         t_range(:,i) = cell2mat(t_tmp);
%     end
%     [long_dot_range, lat_dot_range, r_dot_range] = vel_cart2sph(points(1,:), points(2,:), points(3,:), t_range(1,:), t_range(2,:), t_range(3,:));

    scatter3(points(1,:),points(2,:),points(3,:),20, 'ko')
end


%% Define transfomation matrices
M_tauW = {-sin(lat).*cos(long), -sin(lat).*sin(long),  cos(lat);
          -sin(long)          , cos(long)           ,         0;
          -cos(lat).*cos(long), -cos(lat).*sin(long), -sin(lat)};
M_Wtau = transpose(M_tauW);

M_WO = {cos(obj.ENVMT.windDirection_rad), sin(obj.ENVMT.windDirection_rad), 0;
        sin(obj.ENVMT.windDirection_rad), -cos(obj.ENVMT.windDirection_rad), 0;
        0, 0, -1};
M_OW = transpose(M_WO);
%% Calculate heading difference
t_rot = cell(checks, 1);
t_rot_tau = cell(checks, 1);
chi_parallel = cell(checks, 1);
for i = 1:checks
    t_rot{i} = t_W{i};% doRodriguesRotation(p_kite_W, p_C_W, t);
    t_rot_tau{i} = mult_cellMatrix(M_tauW , t_rot{i});
    chi_parallel{i} = atan2(t_rot_tau{i}{2}, t_rot_tau{i}{1});
end
%% Velocity 
M_AbarO = {cos(chi_a) .* cos(gamma_a), sin(chi_a) .* cos(gamma_a), -sin(gamma_a); 
          -sin(chi_a),                 cos(chi_a),                             0; 
           cos(chi_a) .* sin(gamma_a), sin(chi_a) .* sin(gamma_a), cos(gamma_a)};
M_OAbar = transpose_cellMatrix(M_AbarO);

v_a_O = mult_cellMatrix(M_OAbar,{Va;0;0});
v_k_O = add_cellMatrix(v_a_O, mult_cellMatrix(obj.v_w_O,obj.v0)); 
v_kite_W = mult_cellMatrix(M_WO,v_k_O); 

chi_k = atan2( v_k_O{2}, v_k_O{1} ); 
gamma_k = -asin( v_k_O{3}./norm_cellVec(v_k_O) ); 

%% calculateCommandedCourse2
gamma_cmd = 0;
    
e1 = cell(checks, 1);
e2 = cell(checks, 1);
e3 = cell(checks, 1);
M_CW = cell(checks, 1);
pos_C = cell(checks, 1);
Delta_chi = cell(checks, 1);
Chi_cmd = cell(checks, 1);
v_cmd_tau = cell(checks, 1);
v_cmd_W = cell(checks, 1);
chi_k_ideal = cell(checks, 1);
gamma_k_ideal = cell(checks, 1);
bearing_O = cell(checks, 1);
for i=1:checks
    e1{i} = element_div_cellMatrix(t_rot{i},norm_cellVec(t_rot{i}));
    e3{i} = element_div_cellMatrix(p_C_W{i},-norm_cellVec(p_C_W{i}));
    e2{i} = cross_cellVec(e3{i}, e1{i});

    M_CW{i} = {e1{i}{1} e1{i}{2} e1{i}{3}; ...
     e2{i}{1} e2{i}{2} e2{i}{3}; ...
     e3{i}{1} e3{i}{2} e3{i}{3}};

    pos_C{i} = mult_cellMatrix(M_CW{i},sub_cellMatrix(p_kite_W,p_C_W{i}));

    Delta_chi{i} = atan2( -sign(pos_C{i}{2}) .* delta{i}, delta0);

    % CMD in the tau frame
    Chi_cmd{i} = chi_parallel{i} + Delta_chi{i};
    
    if Chi_cmd{i} > pi
        Chi_cmd{i} = -pi + mod( Chi_cmd{i}, pi );
    elseif Chi_cmd{i} < -pi
        Chi_cmd{i} = pi + mod(Chi_cmd{i}, -pi);
    end
    
    v_cmd_tau{i} = mult_cellMatrix({cos(Chi_cmd{i})*cos(gamma_cmd); 
         sin(Chi_cmd{i})*cos(gamma_cmd); 
        -sin(gamma_cmd)}, {Va});

    v_cmd_W{i} = mult_cellMatrix(M_Wtau,v_cmd_tau{i});
    bearing_O{i} =  mult_cellMatrix(M_OW, v_cmd_W{i});
    chi_k_ideal{i} = atan2( bearing_O{i}{2}, bearing_O{i}{1} ); 
    gamma_k_ideal{i} = -asin( bearing_O{i}{3}./norm_cellVec(bearing_O{i}) ); 

%     chi_cmd_W{i} = atan2( v_cmd_W{i}{2}, v_cmd_W{i}{1});
%     gamma_cmd_W{i} = -asin( v_cmd_W{i}{3}./norm_cellVec(v_cmd_W{i}) ); 
end
% chi_W = atan2( v_kite_W{2}, v_kite_W{1});
% gamma_W = -asin( v_kite_W{3}./norm_cellVec(v_kite_W) ); 

%%
% t_norm = norm_cellVec(t);
% v_kite_W_norm = norm_cellVec(v_kite_W);

% Check which direction leads to a better distance
% distance1 = norm_cellVec(sub_cellMatrix(element_div_cellMatrix(t,t_norm), element_div_cellMatrix(v_kite_W, v_kite_W_norm)));
% distance2 = norm_cellVec(sub_cellMatrix(element_div_cellMatrix(t,-t_norm), element_div_cellMatrix(v_kite_W, v_kite_W_norm)));
chi_delta = cell(checks, 1);
gamma_delta = cell(checks, 1);
distance = cell(checks, 1);
for i=1:checks
    chi_delta{i} = abs(chi_k_ideal{i} - chi_k);
    chi_delta{i} = (chi_delta{i} > pi) .* (2*pi- chi_delta{i}) + (chi_delta{i} <= pi) .* chi_delta{i};
    gamma_delta{i} = abs(gamma_k_ideal{i} - gamma_k);
    gamma_delta{i} = (gamma_delta{i} > pi) .* (2*pi- gamma_delta{i}) + (gamma_delta{i} <= pi) .* gamma_delta{i};
    
    if length(delta{i}) > 1 && normalize
        delta{i} = delta{i}/(1e-99+max(abs(delta{i}(:))));
    end
    if ndims(delta{i}) < 2
        extraOuts.distance_vec{i} = [delta{i}' chi_delta{i}'/(2*pi) gamma_delta{i}'/(2*pi)];
    end

    if isfield(extraArgs, 'distanceOnly') && extraArgs.distanceOnly
        distance{i} = delta{i}(:)';
    elseif isfield(extraArgs, 'headingOnly') && extraArgs.headingOnly
        distance{i} = max(chi_delta{i}/(2*pi), gamma_delta{i}/(2*pi));
        %distance{i} = chi_delta{i};
        distance{i} = distance{i}(:)';
    else
        distance{i} = max(delta{i}, max(chi_delta{i}/(2*pi), gamma_delta{i}/(2*pi)));
        distance{i} = distance{i}(:)';
    end
end

[distance, I] = min(cell2mat(distance), [], 1);
if length(state{1}) > length(state{5})
    distance = reshape(distance, size(delta{1}));
else
    distance = reshape(distance, size(chi_delta{1}));
end

if length(I) == 1
    v_a_cmd_O = sub_cellMatrix(bearing_O{I}, mult_cellMatrix(obj.v_w_O,obj.v0));
    extraOuts.chi_cmd_a = atan2( v_a_cmd_O{2}, v_a_cmd_O{1});
    extraOuts.gamma_cmd_a = -asin( v_a_cmd_O{3}./norm_cellVec(v_a_cmd_O)); 
elseif checks == 1
    v_a_cmd_O = sub_cellMatrix(bearing_O{1}, mult_cellMatrix(obj.v_w_O,obj.v0));
    extraOuts.chi_cmd_a = atan2( v_a_cmd_O{2}, v_a_cmd_O{1});
    extraOuts.gamma_cmd_a = -asin( v_a_cmd_O{3}./norm_cellVec(v_a_cmd_O)); 
else
    extraOuts = [];
end

sol_final = 0;
p_C_W_final = {0;0;0};
v_cmd_W_final = {0;0;0};
t_W_final = {0;0;0};
for i=1:checks
    if any(I(nonconverged_points{i})==i)
        warning('solution point not found')
    end
    sol_final = sol_final + sol{i}(:) .* (I==i)';

    p_C_W_final{1} = p_C_W_final{1} + p_C_W{i}{1}(:) .* (I==i)';
    p_C_W_final{2} = p_C_W_final{2} + p_C_W{i}{2}(:) .* (I==i)';
    p_C_W_final{3} = p_C_W_final{3} + p_C_W{i}{3}(:) .* (I==i)';

    v_cmd_W_final{1} = v_cmd_W_final{1} + v_cmd_W{i}{1}(:) .* (I==i)';
    v_cmd_W_final{2} = v_cmd_W_final{2} + v_cmd_W{i}{2}(:) .* (I==i)';
    v_cmd_W_final{3} = v_cmd_W_final{3} + v_cmd_W{i}{3}(:) .* (I==i)';
    
    t_W_final{1} = t_W_final{1} + t_W{i}{1}(:) .* (I==i)';
    t_W_final{2} = t_W_final{2} + t_W{i}{2}(:) .* (I==i)';
    t_W_final{3} = t_W_final{3} + t_W{i}{3}(:) .* (I==i)';
end
sol = reshape(sol_final, size(distance));

p_C_W = cell(3,1);
p_C_W{1} = reshape(p_C_W_final{1}, size(distance));
p_C_W{2} = reshape(p_C_W_final{2}, size(distance));
p_C_W{3} = reshape(p_C_W_final{3}, size(distance));

v_cmd_W = cell(3,1);
v_cmd_W{1} = reshape(v_cmd_W_final{1}, size(distance));
v_cmd_W{2} = reshape(v_cmd_W_final{2}, size(distance));
v_cmd_W{3} = reshape(v_cmd_W_final{3}, size(distance));

t_W = cell(3,1);
t_W{1} = reshape(t_W_final{1}, size(distance));
t_W{2} = reshape(t_W_final{2}, size(distance));
t_W{3} = reshape(t_W_final{3}, size(distance));
if visualize_distance
    plot3(p_C_W{1}, p_C_W{2}, p_C_W{3}, 'kh')
    if (isfield(extraArgs, 'LongLatState') && extraArgs.LongLatState) || obj.LongLatState
        plot3(p_kite_W{1}, p_kite_W{2}, p_kite_W{3}, 'gp')
    else
        plot3(p_kite_W{1}, p_kite_W{2}, p_kite_W{3}, 'rp')
    end
    
    if isfield(extraArgs, 'visualizeVec') && ~extraArgs.visualizeVec
        return
    end
    
    % translate
    quiver3(p_kite_W{1},p_kite_W{2},p_kite_W{3},v_cmd_W{1},v_cmd_W{2},v_cmd_W{3})

    % plot current velocity
    if length(p_kite_W{1}) == length(v_kite_W{1})
        quiver3(p_kite_W{1},p_kite_W{2},p_kite_W{3},v_kite_W{1},v_kite_W{2},v_kite_W{3})
    end
    
    % visualize velocity at tangent of desired point
    quiver3(p_C_W{1},p_C_W{2},p_C_W{3},t_W{1},t_W{2},t_W{3})
end
end
