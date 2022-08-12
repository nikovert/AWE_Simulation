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

function [dx, F_rest, max_F_tether] = dynamics(obj, ~, x, u, d)
% function dx = dynamics(t, x, u)
%     Dynamics of the 3 DOF AWE system
%
% STATES
%     long    - longitude defined in the wind reference frame 
%     lat     - latitude defined in the wind reference frame
%     h_tau   - radial distance to the origin in the wind reference frame
%     va      - apparent wind speed in the A_bar (rotated aerodynamic frame) frame
%     chi_a   - course angle in the A_bar frame
%     gamma_a - path angle in the A_bar frame
%     Delta t - tether segment difference
%
% INPUTS
%     alpha_a - aerodynamic bank angle defined in the A_bar frame
%     mu_a    - aerodynamic bank angle
%     
%     The sideslip angle (beta_a) is set to zero for the 3 DOF model
%
% obj.ENVMT     - environmental variable
% obj.AIRCRAFT  - Aircarft variables
%

num = false;
if ~iscell(x)
    x = num2cell(x, 2)';
    num = true;
end

if ~iscell(u)
    u = num2cell(u, 2)';
end

if nargin < 5
    d = {0; 0; 0; 0};
end

if ~iscell(d)
    d = num2cell(d, 2)';
end
%% Define states and inputs (remove this step later for speed up)
h_tau    = x{3};
va       = x{4};
chi_a    = x{5} * obj.a0;
gamma_a  = x{6} * obj.a0;

if isempty(obj.curve_direction)
    direction = 1;
else
    direction = obj.curve_direction;
end

if ~isempty(obj.Lem)
    extraArgs.Lem = obj.Lem;
end

extraArgs.direction = direction;

if obj.LongLatState
    long = x{1}* obj.a0;
    lat  = x{2}* obj.a0;
else
    s        = x{1}* obj.a0;
    sigma    = x{2};
    [long, lat, t_tau, t_rot_tau, t_W, t_rot_W] = getLongLat(s, sigma, obj.h0 * h_tau, extraArgs);
end

if length(x) > 6
    tether_diff  = x{7};
else
    tether_diff = 0;
end

alpha    = u{1};
beta     = 0;
mu_a     = u{2};

% Rotational Rates
p = 0; 
q = 0;
r = 0; 

% Hardcoded for now
delta_a = 0; 
delta_e = -0.092174717977942;
delta_r = 0; 

% Replace with disturbance later on
xi = obj.ENVMT.windDirection_rad;
if obj.v0 ~= 1
    warning('wind shear model does not use normalisation')
end
v_w_O = wind_shear(h_tau .* sin(lat) ,xi, obj.base_windspeed);
%% Define rotation matrix
M_AbarA = {1,0,0; 
           0, cos(mu_a), -sin(mu_a); 
           0, sin(mu_a), cos(mu_a)};
% M_AAbar = transpose_cellMatrix(M_AbarA);

M_AB = {cos(alpha) .* cos(beta), sin(beta), sin(alpha) .* cos(beta); 
    -cos(alpha) .* sin(beta), cos(beta), -sin(alpha) .* sin(beta); 
    -sin(alpha), 0, cos(alpha)};
%% Force coefficients 
% Cx calculation
Cx_0 = obj.AIRCRAFT.Cx_0_0 + obj.AIRCRAFT.Cx_0_alpha * alpha + obj.AIRCRAFT.Cx_0_alpha2 * alpha.^2; 
%Cx_q = obj.AIRCRAFT.Cx_q_0 + obj.AIRCRAFT.Cx_q_alpha * alpha + obj.AIRCRAFT.Cx_q_alpha2 * alpha.^2; 
Cx_deltaE  = obj.AIRCRAFT.Cx_deltaE_0  + obj.AIRCRAFT.Cx_deltaE_alpha * alpha + obj.AIRCRAFT.Cx_deltaE_alpha2 * alpha.^2; 
%Cx = Cx_0 + Cx_q*obj.AIRCRAFT.c*q./(2*va*obj.v0) + Cx_deltaE*delta_e; 
Cx = Cx_0 + Cx_deltaE*delta_e; 

% Cy calculation
%Cy_beta = obj.AIRCRAFT.Cy_beta_0 + obj.AIRCRAFT.Cy_beta_alpha * alpha + obj.AIRCRAFT.Cy_beta_alpha2 * alpha.^2; 
%Cy_p = obj.AIRCRAFT.Cy_p_0  + obj.AIRCRAFT.Cy_p_alpha * alpha + obj.AIRCRAFT.Cy_p_alpha2 * alpha.^2;  
%Cy_r = obj.AIRCRAFT.Cy_r_0  + obj.AIRCRAFT.Cy_r_alpha * alpha + obj.AIRCRAFT.Cy_r_alpha2 * alpha.^2;  
Cy_deltaA =obj.AIRCRAFT.Cy_deltaA_0 + obj.AIRCRAFT.Cy_deltaA_alpha * alpha + obj.AIRCRAFT.Cy_deltaA_alpha2 * alpha.^2;  
Cy_deltaR =obj.AIRCRAFT.Cy_deltaR_0 + obj.AIRCRAFT.Cy_deltaR_alpha * alpha + obj.AIRCRAFT.Cy_deltaR_alpha2 * alpha.^2;  
%Cy = Cy_beta*beta + Cy_p*obj.AIRCRAFT.b*p./(2*va*obj.v0) + Cy_r*obj.AIRCRAFT.b*r./(2*va*obj.v0) + Cy_deltaA*delta_a + Cy_deltaR*delta_r;
Cy = Cy_deltaA*delta_a + Cy_deltaR*delta_r;

% Cz calculation
Cz_0 = obj.AIRCRAFT.Cz_0_0 + obj.AIRCRAFT.Cz_0_alpha * alpha + obj.AIRCRAFT.Cz_0_alpha2 * alpha.^2; 
%Cz_q = obj.AIRCRAFT.Cz_q_0 + obj.AIRCRAFT.Cz_q_alpha * alpha + obj.AIRCRAFT.Cz_q_alpha2 * alpha.^2; 
Cz_deltaE  = obj.AIRCRAFT.Cz_deltaE_0 + obj.AIRCRAFT.Cz_deltaE_alpha * alpha + obj.AIRCRAFT.Cz_deltaE_alpha2 * alpha.^2; 
%Cz = Cz_0 + Cz_q*obj.AIRCRAFT.c*q./(2*va*obj.v0) + Cz_deltaE*delta_e; 
Cz = Cz_0 + Cz_deltaE*delta_e; 
%% Translation Dynamics 
% if num || (isempty(obj.s_dot) ||isempty(obj.sigma_dot) ||isempty(obj.long_dot) || isempty(obj.lat_dot) || isempty(obj.h_tau_dot))
    
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
    
    obj.long_dot = (v_k_tau{2})./( abs(h_tau).*cos(lat));
    obj.lat_dot = (v_k_tau{1})./abs(h_tau);
    obj.h_tau_dot = -(v_k_tau{3});
    if ~obj.LongLatState
        Lem = extraArgs.Lem;
        Lem.a = Lem.a./h_tau;
        Lem.b = Lem.b./h_tau;
        
%         dlongds = @(s) ( Lem.b.^3 .* cos(s).*(2*Lem.a.^2-Lem.a.^2 .* cos(s).^2+Lem.b.^2)./(Lem.a.^2 .* cos(s).^2+Lem.b.^2).^2 );
%         dlatds  = @(s) ((cos(s).^2 .* (Lem.a.^3 .* Lem.b.^2+2*Lem.a .* Lem.b.^4) - Lem.a .* Lem.b.^4)./(Lem.a.^2 .* cos(s).^2+Lem.b.^2).^2 );
% 
%         dxds = @(s) h_tau * (cos(lat) .* dlongds(s) - sin(lat) .* dlatds(s) .* long);
%         dyds = @(s) h_tau * dlatds(s);
%         s_dot = ([dxds(s), dyds(s)] * [v_k_tau{1}; v_k_tau{2}])/ (norm([dxds(s), dyds(s)]) * h_tau);
%         
        obj.s_dot     = cell2mat(mult_cellMatrix(transpose(t_W), mult_cellMatrix(v_k_W ,obj.v0))) ...
            ./ (direction * norm_cellVec(t_W));
        
        obj.s_dot = 2*pi/945 * obj.s_dot;
        
        obj.sigma_dot = cell2mat(mult_cellMatrix(transpose(t_rot_W), mult_cellMatrix(v_k_W ,obj.v0))) ...
            ./ (norm_cellVec(t_rot_W));
    end
% end
if ~obj.LongLatState
    s_dot     = obj.s_dot;
    sigma_dot = obj.sigma_dot;
end
% if sigma_dot > 0
%    warning('bad sigma') 
% end
long_dot  = obj.long_dot;
lat_dot   = obj.lat_dot;
h_tau_dot = obj.h_tau_dot;


%% For checks
% if ~obj.LongLatState
%     Lem.a = obj.Lem.a ./ (obj.h0*h_tau);
%     Lem.b = obj.Lem.b ./ (obj.h0*h_tau);
%     Lem.phi0 = obj.Lem.phi0;
%     
%     M_WP = [cos(Lem.phi0),0, -sin(Lem.phi0);0, 1, 0; sin(Lem.phi0),0, cos(Lem.phi0)];
%     
%     [t,DtDs,L_P, dLds_P] = getBoothInfos2(s,Lem, direction);
%     
%     [vx_P, vy_P, vz_P] = vel_sph2cart(L_P{1}, L_P{2}, h_tau, dLds_P{1}, dLds_P{2}, h_tau_dot);
%     V = M_WP * [vx_P, vy_P, vz_P]';
%     
%     [pos_W_x,pos_W_y,pos_W_z] = sph2cart(long,lat,obj.h0*h_tau);
%     [dLds_P_long, dLds_P_lat, r_dot] = vel_cart2sph(pos_W_x, pos_W_y, pos_W_z, V(1), V(2), V(3));
%      
%     long_dot_alt = dLds_P_long * s_dot;
%     lat_dot_alt  = dLds_P_lat * s_dot;
% 
% %     s_dot = (long_dot/dLds_P_long + lat_dot/dLds_P_lat)/2;
% %     if abs(long_dot/dLds_P_long - lat_dot/dLds_P_lat) > 0.01
% %         warning('missmatch')
% %     end
% end
%% Calculate Tether forces
% if num || (isempty(obj.pos_W))
%     [pos_W_x,pos_W_y,pos_W_z] = sph2cart(long,lat,obj.h0*h_tau);
%     pos_W{1,1} = pos_W_x;
%     pos_W{2,1} = pos_W_y;
%     pos_W{3,1} = pos_W_z;
%     obj.pos_W = pos_W;
% else
%     pos_W = obj.pos_W;
% end
%% Calculate Aerodynamic Forces and Moments
F_a_B = mult_cellMatrix({0.5 * obj.AIRCRAFT.S_ref * 1.225*(va*obj.v0).^2}, {Cx;Cy; Cz});
F_a_A = mult_cellMatrix(M_AB, F_a_B); % checked mon 16 Aug

%% All Forces except Aero
% if num || isempty(obj.F_rest)

%     if false % simple tether
%         pos_W_norm = norm_cellVec(pos_W);
%         F_t_W = mult_cellMatrix(element_div_cellMatrix(pos_W,{-pos_W_norm}),obj.Ft_set);
%         l_s_dot = h_tau_dot./(obj.T.n_t_p+1);
%         max_F_tether = 0;
%     else
        x_W{1} = long;
        x_W{2} = lat;
        x_W{3} = h_tau*obj.h0;
        x_W{4} = long_dot * obj.v0/obj.h0;
        x_W{5} = lat_dot * obj.v0/obj.h0;
        x_W{6} = h_tau_dot*obj.v0;
        x_W{7} = tether_diff;
        if nargout > 2
            [f_kite_G, ~, max_F_tether] = obj.tether_forces(x_W);
        else
            [f_kite_G, ~] = obj.tether_forces(x_W);
        end
        F_t_W = mult_cellMatrix({-1}, f_kite_G);
%     end
    
%     F_t_W_norm = norm_cellVec(F_t_W);
%     F_t_W = mult_cellMatrix(element_div_cellMatrix(F_t_W, 1e-99+F_t_W_norm), {min(F_t_W_norm, obj.F_T_max*1e+03)});
%     
    f_kite_G_O = mult_cellMatrix(M_OW,F_t_W);
    if obj.skipTether
        F_rest = mult_cellMatrix(M_AbarO, {0;0;obj.AIRCRAFT.mass * obj.ENVMT.g});
    else
        F_rest = mult_cellMatrix(M_AbarO, add_cellMatrix({0;0;obj.AIRCRAFT.mass * obj.ENVMT.g},f_kite_G_O));
    end
    obj.F_rest = F_rest;
    %obj.final_seg_diff_dot = final_seg_diff_dot;
% else
%     max_F_tether = 0;
%     F_rest = obj.F_rest;
% end
%tether_diff_dot = cell2mat(obj.final_seg_diff_dot) - l_s_dot;
tether_diff_dot = d{1};
if obj.ignoreTetherDiff
    tether_diff_dot = 0;
end
%% Tether Drag
% va_tau = mult_cellMatrix(M_tauW, mult_cellMatrix(M_WO,v_a_O));
% Fd_tau = mult_cellMatrix({-1/8*1.225*obj.T.CD_tether*pos_W_norm*T.d_tether* sqrt( va_tau{1}.^2 + va_tau{2}.^2 )} , {Va_tau{1}; Va_tau{2}; 0});
% 
% Fd_W = M_tauW'*Fd_tau; %transformFromTautoW( long, lat, Fd_tau);
% Fd_O = transformFromWtoO(  windDirection_rad, Fd_W);
% Ftd_B2 = M_OB'*Fd_O;% transformFromAbar2B(mu_a, alpha, beta, Fd_Abar );
% 
% Ftd_A2 = M_AB*Ftd_B2;
%else
%simplified
% if num || (isempty(obj.F_tether_drag_Abar))
%     pos_W_norm = norm_cellVec(pos_W);
%     F_tether_drag_Abar = mult_cellMatrix({-1/8*1.225*obj.T.CD_tether*pos_W_norm*obj.T.d_tether.*va*obj.v0},{va*obj.v0;0;0});
%     obj.F_tether_drag_Abar = F_tether_drag_Abar;
% else
%     F_tether_drag_Abar = obj.F_tether_drag_Abar;
% end
%% Combine all forces
% Ftot_Abar = add_cellMatrix(mult_cellMatrix(M_AbarA,add_cellMatrix(F_a_A,F_rest)),F_tether_drag_Abar);
Ftot_Abar = add_cellMatrix(mult_cellMatrix(M_AbarA,F_a_A),F_rest);
%Ftot_Abar = add_cellMatrix(add_cellMatrix(mult_cellMatrix(M_AbarA,F_a_A),Fg_Abar), mult_cellMatrix(M_AbarO,kite_tether_force(:, end)));

%%
va_dot = 1/obj.AIRCRAFT.mass * Ftot_Abar{1} * obj.h0/obj.v0^2;
chi_a_dot  = 1/obj.AIRCRAFT.mass * 1./(va*obj.v0.*cos(gamma_a)) .* Ftot_Abar{2} * obj.h0/obj.v0;
gamma_a_dot = 1/obj.AIRCRAFT.mass * -1./(va*obj.v0) .* Ftot_Abar{3} * obj.h0/obj.v0;

%%
if num
    dx = nan(obj.nx,length(h_tau));
    if obj.LongLatState
        dx(1,:) = long_dot * obj.v0/obj.h0/obj.a0;
        dx(2,:) = lat_dot * obj.v0/obj.h0/obj.a0;
    else
        dx(1,:) = s_dot* obj.v0/obj.h0/obj.a0;
        dx(2,:) = sigma_dot* obj.v0/obj.h0;
    end
    dx(3,:) = h_tau_dot;
    dx(4,:) = va_dot;
    dx(5,:) = chi_a_dot/obj.a0;
    dx(6,:) = gamma_a_dot/obj.a0;
    if(obj.nx > 6)
        dx(7,:) = tether_diff_dot * obj.v0/obj.h0;
    end
    if ~isreal(dx)
        warning('complex dx')
    end
else 
    dx = cell(obj.nx,1);
    if obj.LongLatState
        dx{1} = long_dot * obj.v0/obj.h0/obj.a0;
        dx{2} = lat_dot * obj.v0/obj.h0/obj.a0;
    else
        dx{1} = s_dot* obj.v0/obj.h0/obj.a0;
        dx{2} = sigma_dot* obj.v0/obj.h0;
    end
    dx{3} = h_tau_dot;
    dx{4} = va_dot;
    dx{5} = chi_a_dot/obj.a0;
    dx{6} = gamma_a_dot/obj.a0;
    if(obj.nx > 6)
        dx{7} = tether_diff_dot* obj.v0/obj.h0;
    end
end

end
