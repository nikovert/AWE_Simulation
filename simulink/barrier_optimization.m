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

function [mu_a, alpha_a] = barrier_optimization(Va, u_ref, A, b, P_AP2, constr)
%BARRIER_OPTIMIZATION compute the optimal compromise betweeen u_safety and
%u_NDI based on a control barrier function
    
    problem.objective = @(u) (u - u_ref)'*(u - u_ref);
    problem.solver = 'fmincon';
    problem.x0 = u_ref;
    problem.lb = [constr.mu_a_min constr.alpha_a_min];
    problem.ub = [constr.mu_a_max constr.alpha_a_max];
    problem.Aineq = [];
    problem.bineq = [];
    problem.Aeq = [];
    problem.beq = [];
    problem.nonlcon = @(u) nonlconstraint(u, Va, A, b, P_AP2);
    problem.options = optimoptions('fmincon', 'Display', 'off','Algorithm','sqp');

    u_robust = fmincon(problem);
    mu_a    =  u_robust(1);
    alpha_a =  u_robust(2);
end

function [c,ceq] = nonlconstraint(u, Va, A, b, P_AP2)
    mu_a    =  u(1);
    alpha =  u(2);
    F_a_Abar = computeFa(Va, alpha, mu_a, P_AP2);

    c = double(-(A'*F_a_Abar+b));
    ceq = 0;
end

function F_a_Abar = computeFa(Va, alpha, mu_a, P_AP2)
    % Rotational Rates
    p =0; 
    q = 0; %-0.1202; 
    r = 0; 
    
    beta = 0; 
    
    delta_a = 0; 
    delta_e = -0.092174717977942; % apparently has quite an impact on the Cx coefficient, hardcoded for now
    delta_r = 0; 
    
    %% Force coefficients 
    % Cx calculation
    Cx_0 = P_AP2.Cx_0_0 + P_AP2.Cx_0_alpha * alpha + P_AP2.Cx_0_alpha2 * alpha^2; 
    Cx_q = P_AP2.Cx_q_0 + P_AP2.Cx_q_alpha * alpha + P_AP2.Cx_q_alpha2 * alpha^2; 
    Cx_deltaE  = P_AP2.Cx_deltaE_0  + P_AP2.Cx_deltaE_alpha * alpha + P_AP2.Cx_deltaE_alpha2 * alpha^2; 
    Cx = Cx_0 + Cx_q*P_AP2.c*q/(2*Va) + Cx_deltaE*delta_e; 
    
    % Cy calculation
    Cy_beta = P_AP2.Cy_beta_0 + P_AP2.Cy_beta_alpha * alpha + P_AP2.Cy_beta_alpha2 * alpha^2; 
    Cy_p = P_AP2.Cy_p_0  + P_AP2.Cy_p_alpha * alpha + P_AP2.Cy_p_alpha2 * alpha^2;  
    Cy_r = P_AP2.Cy_r_0  + P_AP2.Cy_r_alpha * alpha + P_AP2.Cy_r_alpha2 * alpha^2;  
    Cy_deltaA =P_AP2.Cy_deltaA_0 + P_AP2.Cy_deltaA_alpha * alpha + P_AP2.Cy_deltaA_alpha2 * alpha^2;  
    Cy_deltaR =P_AP2.Cy_deltaR_0 + P_AP2.Cy_deltaR_alpha * alpha + P_AP2.Cy_deltaR_alpha2 * alpha^2;  
    Cy = Cy_beta*beta + Cy_p*P_AP2.b*p/(2*Va) + Cy_r*P_AP2.b*r/(2*Va) + Cy_deltaA*delta_a + Cy_deltaR*delta_r;
    
    % Cz calculation
    Cz_0 = P_AP2.Cz_0_0 + P_AP2.Cz_0_alpha * alpha + P_AP2.Cz_0_alpha2 * alpha.^2; 
    Cz_q = P_AP2.Cz_q_0 + P_AP2.Cz_q_alpha * alpha + P_AP2.Cz_q_alpha2 * alpha^2; 
    Cz_deltaE  = P_AP2.Cz_deltaE_0 + P_AP2.Cz_deltaE_alpha * alpha + P_AP2.Cz_deltaE_alpha2 * alpha^2; 
    Cz = Cz_0 + Cz_q*P_AP2.c*q/(2*Va) + Cz_deltaE*delta_e; 
    
    F_a_B = 0.5*Va^2*P_AP2.S_ref*1.225* [Cx;Cy; Cz];
    
    F_a_Abar = transformFromBtoAbar( mu_a, alpha, beta,  F_a_B );
end