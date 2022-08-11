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
% :Revision: 04-December-2022
% :Author: Nikolaus Vertovec (nikolaus.vertovec@eng.ox.ac.uk)

function [mu_a, alpha_a, g_opt] = barrier_optimization(s, sigma, h_tau, Va, chi_a, gamma_a, tether_diff, d, F_t_W, u_safety, u_NDI, T, ENVMT, P_AP2, base_windspeed, params, PHI_BOOTH, Ft_max, g_opt_prev)
%BARRIER_OPTIMIZATION compute the optimal compromise betweeen u_safety and
%u_NDI based on a control barrier function
    
    x0 = double([s, sigma, h_tau, Va, chi_a, gamma_a, tether_diff, F_t_W]');
    %dx = dynamics(0, x0, u, d, PHI_BOOTH, T, ENVMT, P_AP2, base_windspeed, params);
    tspan =double([0 0.1]);
    options = odeset('RelTol',1e-3,'AbsTol',1e-4, 'MaxStep', 0.1);

    problem.objective = @(g) g;
    problem.solver = 'fmincon';
    problem.x0 = g_opt_prev;
    problem.lb = 0;
    problem.ub = 1;
    problem.Aineq = [];
    problem.bineq = [];
    problem.Aeq = [];
    problem.beq = [];
    problem.nonlcon = @(g) nonlconstraint(g, u_safety, u_NDI, d, x0, tspan, Ft_max, PHI_BOOTH, T, ENVMT, P_AP2, base_windspeed, params, options);
    problem.options = optimoptions('fmincon', 'Display', 'off');

    g_opt = fmincon(problem);
    mu_a    =  (1-g_opt)*u_NDI(1) + g_opt*u_safety(1);
    alpha_a =  (1-g_opt)*u_NDI(2) + g_opt*u_safety(2);
end

function dx  = dynamics(t, x, u, d, PHI_BOOTH, T, ENVMT, P_AP2, base_windspeed, params)
    
    s           = x(1);
    sigma       = x(2);
    h_tau       = x(3);
    va          = x(4);
    chi_a       = x(5);
    gamma_a     = x(6);
    tether_diff = x(7);
    F_t_W       = x(8);
    % 
    alpha       = u(1);
    mu_a        = u(2);
    xi          = ENVMT.windDirection_rad;

    direction = params.direction;
    [long, lat, t_tau, t_rot_tau] = getLongLat(s, sigma, h_tau, PHI_BOOTH, params);
    
    vw = base_windspeed * log(3.281 * h_tau .* sin(lat)/0.15)/log(20/0.15);
    v_w_O = [cos(xi) * vw; sin(xi) * vw; 0];

    beta     = 0;
    % Rotational Rates
    p = 0; 
    q = 0;
    r = 0; 
    
    % Hardcoded for now
    delta_a = 0; 
    delta_e = -0.092174717977942;
    delta_r = 0; 
    
    %% Define rotation matrix
    M_AbarA = [1,0,0; 
               0, cos(mu_a), -sin(mu_a); 
               0, sin(mu_a), cos(mu_a)];
    
    M_AB = [cos(alpha) .* cos(beta), sin(beta), sin(alpha) .* cos(beta); 
        -cos(alpha) .* sin(beta), cos(beta), -sin(alpha) .* sin(beta); 
        -sin(alpha), 0, cos(alpha)];
    
    %% Force coefficients 
    % Cx calculation
    Cx_0 = P_AP2.Cx_0_0 + P_AP2.Cx_0_alpha * alpha + P_AP2.Cx_0_alpha2 * alpha.^2; 
    Cx_q = P_AP2.Cx_q_0 + P_AP2.Cx_q_alpha * alpha + P_AP2.Cx_q_alpha2 * alpha.^2; 
    Cx_deltaE  = P_AP2.Cx_deltaE_0  + P_AP2.Cx_deltaE_alpha * alpha + P_AP2.Cx_deltaE_alpha2 * alpha.^2; 
    Cx = Cx_0 + Cx_q*P_AP2.c*q./(2*va) + Cx_deltaE*delta_e; 
    
    % Cy calculation
    Cy_beta = P_AP2.Cy_beta_0 + P_AP2.Cy_beta_alpha * alpha + P_AP2.Cy_beta_alpha2 * alpha.^2; 
    Cy_p = P_AP2.Cy_p_0  + P_AP2.Cy_p_alpha * alpha + P_AP2.Cy_p_alpha2 * alpha.^2;  
    Cy_r = P_AP2.Cy_r_0  + P_AP2.Cy_r_alpha * alpha + P_AP2.Cy_r_alpha2 * alpha.^2;  
    Cy_deltaA = P_AP2.Cy_deltaA_0 + P_AP2.Cy_deltaA_alpha * alpha + P_AP2.Cy_deltaA_alpha2 * alpha.^2;  
    Cy_deltaR = P_AP2.Cy_deltaR_0 + P_AP2.Cy_deltaR_alpha * alpha + P_AP2.Cy_deltaR_alpha2 * alpha.^2;  
    Cy = Cy_beta*beta + Cy_p*P_AP2.b*p./(2*va) + Cy_r*P_AP2.b*r./(2*va) + Cy_deltaA*delta_a + Cy_deltaR*delta_r;
    
    % Cz calculation
    Cz_0 = P_AP2.Cz_0_0 + P_AP2.Cz_0_alpha * alpha + P_AP2.Cz_0_alpha2 * alpha.^2; 
    Cz_q = P_AP2.Cz_q_0 + P_AP2.Cz_q_alpha * alpha + P_AP2.Cz_q_alpha2 * alpha.^2; 
    Cz_deltaE  = P_AP2.Cz_deltaE_0 + P_AP2.Cz_deltaE_alpha * alpha + P_AP2.Cz_deltaE_alpha2 * alpha.^2; 
    Cz = Cz_0 + Cz_q*P_AP2.c*q./(2*va) + Cz_deltaE*delta_e; 
    
    %% Translation Dynamics 
    M_tauW = [-sin(lat) .* cos(long), -sin(lat) .* sin(long), cos(lat);
              -sin(long)            , cos(long)             , 0;
              -cos(lat) .* cos(long), -cos(lat) .* sin(long), -sin(lat) ];
    
    M_AbarO = [cos(chi_a) .* cos(gamma_a), sin(chi_a) .* cos(gamma_a), -sin(gamma_a); 
              -sin(chi_a),                 cos(chi_a),                             0; 
               cos(chi_a) .* sin(gamma_a), sin(chi_a) .* sin(gamma_a), cos(gamma_a)];
    M_OAbar = M_AbarO';
    
    M_OW = [cos(xi),  sin(xi),  0; 
            sin(xi), -cos(xi),  0; 
            0,              0, -1];
    M_WO = M_OW; 
    
    v_a_O = M_OAbar * [va;0;0];
    v_k_O = v_a_O + v_w_O; 
    v_k_W = M_WO * v_k_O; 
    v_k_W = v_k_W + d(2:4);
    
    v_k_tau = M_tauW * v_k_W;
    
    
    s_dot     = 2*pi/945 * t_tau'*v_k_tau./ (direction*norm(t_tau));
    sigma_dot = t_rot_tau'*v_k_tau./ (norm(t_tau));
    
    long_dot = (v_k_tau(2))./( abs(h_tau).*cos(lat));
    lat_dot = (v_k_tau(1))./abs(h_tau);
    h_tau_dot = -(v_k_tau(3));
    
    [pos_W_x,pos_W_y,pos_W_z] = sph2cart(long, lat, h_tau);
    pos_W = [pos_W_x; pos_W_y; pos_W_z];
    %% Calculate Aerodynamic Forces and Moments
    F_a_B = 0.5 * P_AP2.S_ref * 1.225*(va).^2 * [Cx;Cy; Cz];
    F_a_A = M_AB * F_a_B;
    
    %% All Forces except Aero
    % x_W = [long; lat; h_tau; long_dot; lat_dot; h_tau_dot; tether_diff];
    % [f_kite_G_approx , final_seg_diff_dot] = tether_forces(x_W, P_AP2, ENVMT, T, base_windspeed);
    %tether_diff_dot = final_seg_diff_dot - l_s_dot;
    % F_t_W_approx = -1 * f_kite_G_approx;
    
    tether_diff_dot = d(1);

    f_kite_G_O = M_OW * F_t_W;
    F_rest_Abar = M_AbarO * ([0;0;P_AP2.mass * ENVMT.g] + f_kite_G_O);
    
    %% Combine all forces
    Ftot_Abar = M_AbarA * F_a_A + F_rest_Abar;
    
    va_dot = 1/P_AP2.mass * Ftot_Abar(1);
    chi_a_dot  = 1/P_AP2.mass * 1/(va .* cos(gamma_a)) .* Ftot_Abar(2);
    gamma_a_dot = 1/P_AP2.mass * -1/(va) .* Ftot_Abar(3);
    dx = double([s_dot, sigma_dot, h_tau_dot, va_dot, chi_a_dot, gamma_a_dot, tether_diff_dot, T.c0 * tether_diff_dot]');
end

function [c,ceq] = nonlconstraint(g, u_safety, u_NDI, d, x0, tspan, Ft_max, PHI_BOOTH, T, ENVMT, P_AP2, base_windspeed, params, options)
    u = (1-g)*u_safety + g*u_NDI;
    [~, x] = ode23(@(t,x) dynamics(t, x, u, d, PHI_BOOTH, T, ENVMT, P_AP2, base_windspeed, params), tspan, x0, options);
    rho = Ft_max-x(:,end);
    h = min(rho);
    c = -h;
    ceq = 0;
end