% Copyright (C) 2023  Nikolaus Vertovec
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the iFTmplied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
% :Revision: 31-January-2023
% :Author: Nikolaus Vertovec (nikolaus.vertovec@eng.ox.ac.uk)
% :Adapted from: Sebastian Rapp (s.rapp@tudelft.nl)

function [x0_sim,u0_sim,  act, aeroModel, base_windspeed, constr,...
    ENVMT, P_AP2, simInit, T, winchParameter,params] = initAllSimParams(lat_in)
simInit.TSIM = 1200; 
simInit.Ts_power_conv_check = 0.1; %Power convergence check sample time
simInit.power_conv_threshold = 500; 
simInit.Ts_vis = 5;
params.a_booth = 0.6; % 0.6
params.b_booth = 200; % 200 also fine
params.phi0_booth =  30*pi/180;
base_windspeed = 9; %
params.meas_delay = 0.015; 
params.control_delay = 0.015; 
params.winch_control_delay = 0.015; 
params.winch_meas_delay = 0.015; 
simInit.skipT = 60; 
params.Ft_set_traction = 1400; 
params.Ft_set_retraction = 500; %500; 
params.winch_control_filter_w0 =1*2*pi; 
params.winch_control_filter_w0 = 2*2*pi; 
params.phi0_dot = 0.1; 
params.eject2retraction = 5*pi/180;% 10*pi/180; 
params.Ft_set_init = 1000; 
% these flags define if the aerodynamic model is perturbed or not and if
% the tether force acts in the CG or not (o.w. it induces a moment)
simInit.perturbed_aero_flag = 0; 
simInit.perturbed_ft_moment_flag = 0; 
simInit.complex_tether_flag = 1; 

params.vr_mean_min = 2; 
params.Ft_set_traction_low = 800; 

u0_sim = [-0.0023   -0.0922   -0.0081 1800]';
x0_sim = [31.3466         0    0.1222    0.1466    0.0429   -1.8953    0.0176   -0.0967    0.0148   -0.0000    1.3963  250.0000]';

x0_sim(11) = lat_in;
x0_sim(12) = 250;
pos_init_W = [cos( x0_sim(10) ) * cos( lat_in );
    sin( x0_sim(10) ) * cos( lat_in );
    sin( lat_in )]*x0_sim(12);
%% Environment struct
ENVMT.g= 9.8066;
ENVMT.rhos= 1.2250;
ENVMT.windDirection_rad= 3.1416;

simInit.pos_O_init = transformFromWtoO(ENVMT.windDirection_rad, pos_init_W)'; %[-200,0,-250];
simInit.tether_inital_lenght = norm( simInit.pos_O_init)+1;
simInit.pos_W_init = transformFromOtoW(ENVMT.windDirection_rad, simInit.pos_O_init');
simInit.long_init = atan2( simInit.pos_W_init(2), simInit.pos_W_init(1) );
simInit.lat_init = asin( simInit.pos_W_init(3)/norm(simInit.pos_W_init) );
simInit.dt= 2.5000e-04;

%% Actuators
% Aileron
act.aileron.w0 = 35; %35 is more realistic
act.aileron.rl = 2;%300*pi/180; %
act.aileron.max_deflection = 20*pi/180;
act.aileron.min_deflection = -20*pi/180;
act.aileron.bias = 0;

% Elevator
act.elevator.w0 = 35; %
act.elevator.rl = 2;%300*pi/180; %
act.elevator.max_deflection = 20*pi/180;
act.elevator.min_deflection = -20*pi/180;
act.elevator.bias = 0;

% Rudder
act.rudder.w0 = 35; %
act.rudder.rl = 2;%300*pi/180; %
act.rudder.max_deflection = 30*pi/180;
act.rudder.min_deflection = -30*pi/180;
act.rudder.bias = 0;

%% Constrainst
constr.mu_a_max= 1.0472;
constr.mu_a_min= -1.0472;
constr.alpha_a_max= 10*pi/180;
constr.alpha_a_min= -0.1047;
constr.max_sideslip = 20*pi/180; 

constr.mu_a_max_retraction= 30*pi/180;% 0.7854;
constr.mu_a_min_retraction= -30*pi/180;%-0.7854;
constr.max_lift= 2000;
constr.max_delta_mu= 0.1745;
constr.mu_k_max= 1.0472;
constr.phi_tau_max= 1.3963;

constr.winchParameter.a_max= 5;
constr.winchParameter.a_min= -5;
constr.winchParameter.v_max= 20;
constr.winchParameter.v_min=-15;
constr.cte_max = 60;
constr.gs_error = 60; 

constr.gamma_retraction_min = -20*pi/180; 
constr.gamma_retraction_max = 20*pi/180; 

constr.max_CL= 2;
constr.min_CL= -2;

constr.F_t_max = 1870; 

constr.vamin = 25; 

constr.max_rates_B = 1; %0.8727;

%% Aircraft parameters
P_AP2.AP_Ts= 0.0100;
P_AP2.delay_g2k= 0.1000;
P_AP2.delay_k2g= 0.1000;
P_AP2.S_ref= 3;
P_AP2.S_wing= 3;
P_AP2.b= 5.5000;
P_AP2.c= 0.5500;
P_AP2.mass= 36.8000;
P_AP2.J=  [25  ,       0 ,  -0.47;
    0 ,  32   ,      0;
    -0.47   ,      0 ,  56];

P_AP2.Jx= 25;
P_AP2.Jy= 32;
P_AP2.Jz= 56;
P_AP2.Jxz= 0.4700;
% P_AP2.Jinv= [0.0400   ,      0  ,  0.0003;
%     0 ,   0.0313      ,   0;
%     0.0003    ,     0 ,   0.0179];
P_AP2.Jinv = inv( P_AP2.J ); 
P_AP2.d_tether= 0.0025;
P_AP2.rho_tether= 0.0046;
P_AP2.Cd_tether= 1.2000;
P_AP2.rho_air= 1.2250;
P_AP2.Cx_0_0= -0.0293;
P_AP2.Cx_0_alpha= 0.4784;
P_AP2.Cx_0_alpha2= 2.5549;
P_AP2.Cx_q_0= -0.6029;
P_AP2.Cx_q_alpha= 4.4124;
P_AP2.Cx_q_alpha2= 0;
P_AP2.Cx_deltaE_0= -0.0106;
P_AP2.Cx_deltaE_alpha= 0.1115;
P_AP2.Cx_deltaE_alpha2= 0;
P_AP2.Cy_beta_0= -0.1855;
P_AP2.Cy_beta_alpha= -0.0299;
P_AP2.Cy_beta_alpha2= 0.0936;
P_AP2.Cy_p_0= -0.1022;
P_AP2.Cy_p_alpha= -0.0140;
P_AP2.Cy_p_alpha2= 0.0496;
P_AP2.Cy_r_0= 0.1694;
P_AP2.Cy_r_alpha= 0.1368;
P_AP2.Cy_r_alpha2= 0;
P_AP2.Cy_deltaA_0= -0.0514;
P_AP2.Cy_deltaA_alpha= -0.0024;
P_AP2.Cy_deltaA_alpha2= 0.0579;
P_AP2.Cy_deltaR_0= 0.1032;
P_AP2.Cy_deltaR_alpha= 0.0268;
P_AP2.Cy_deltaR_alpha2= -0.1036;
P_AP2.Cz_0_0= -0.5526;
P_AP2.Cz_0_alpha= -5.0676;
P_AP2.Cz_0_alpha2= 5.7736;
P_AP2.Cz_q_0= -7.5560;
P_AP2.Cz_q_alpha= 0.1251;
P_AP2.Cz_q_alpha2= 6.1486;
P_AP2.Cz_deltaE_0= -0.3150;
P_AP2.Cz_deltaE_alpha= -0.0013;
P_AP2.Cz_deltaE_alpha2= 0.2923;
P_AP2.Cm_0_0= -0.0307;
P_AP2.Cm_0_alpha= -0.6027;
P_AP2.Cm_0_alpha2= 0;
P_AP2.Cm_q_0= -11.3022;
P_AP2.Cm_q_alpha= -0.0026;
P_AP2.Cm_q_alpha2= 5.2885;
P_AP2.Cm_deltaE_0= -1.0427;
P_AP2.Cm_deltaE_alpha= -0.0061;
P_AP2.Cm_deltaE_alpha2= 0.9974;
P_AP2.Cl_beta_0= -0.0630;
P_AP2.Cl_beta_alpha= -3.0000e-04;
P_AP2.Cl_beta_alpha2= 0.0312;
P_AP2.Cl_p_0= -0.5632;
P_AP2.Cl_p_alpha= -0.0247;
P_AP2.Cl_p_alpha2= 0.2813;
P_AP2.Cl_r_0= 0.1811;
P_AP2.Cl_r_alpha= 0.6448;
P_AP2.Cl_r_alpha2= 0;
P_AP2.Cl_deltaA_0= -0.2489;
P_AP2.Cl_deltaA_alpha= -0.0087;
P_AP2.Cl_deltaA_alpha2= 0.2383;
P_AP2.Cl_deltaR_0= 0.0044;
P_AP2.Cl_deltaR_alpha= -0.0013;
P_AP2.Cl_deltaR_alpha2= 0;
P_AP2.Cn_beta_0= 0.0577;
P_AP2.Cn_beta_alpha= -0.0849;
P_AP2.Cn_beta_alpha2= 0;
P_AP2.Cn_p_0= -0.0565;
P_AP2.Cn_p_alpha= -0.9137;
P_AP2.Cn_p_alpha2= 0;
P_AP2.Cn_r_0= -0.0553;
P_AP2.Cn_r_alpha= 0.0290;
P_AP2.Cn_r_alpha2= 0.0257;
P_AP2.Cn_deltaA_0= 0.0190;
P_AP2.Cn_deltaA_alpha= -0.1147;
P_AP2.Cn_deltaA_alpha2= 0;
P_AP2.Cn_deltaR_0= -0.0404;
P_AP2.Cn_deltaR_alpha= -0.0117;
P_AP2.Cn_deltaR_alpha2= 0.0409;
P_AP2.CL0= 0.4687;
P_AP2.CL_alpha= 4.5619;
P_AP2.g= 9.8100;
P_AP2.sampling_rate_fcs= 100;
P_AP2.k_motor= 80;
P_AP2.S_prop= 0.2027;
P_AP2.C_prop= 1;
P_AP2.maxBandwidth= 35;

%% Aerodynamic model that the controller knows
% Aerodynamic uncertainties
aeroModel.CL0 = 1.*P_AP2.CL0;
aeroModel.CL_alpha = 1.*P_AP2.CL_alpha;

% first is scale second bias
s = 1;
b = 0;
aeroModel.Cx_0_0  = addAdditiveAndMultiplUnc(P_AP2.Cx_0_0,s, b);
aeroModel.Cx_0_alpha   = addAdditiveAndMultiplUnc(P_AP2.Cx_0_alpha,s, b);
aeroModel.Cx_0_alpha2  = addAdditiveAndMultiplUnc(P_AP2.Cx_0_alpha2,s,b);

aeroModel.Cy_beta_0  = addAdditiveAndMultiplUnc(P_AP2.Cy_beta_0,s,b);
aeroModel.Cy_beta_alpha  = addAdditiveAndMultiplUnc(P_AP2.Cy_beta_alpha,s, b);
aeroModel.Cy_beta_alpha2 = addAdditiveAndMultiplUnc(P_AP2.Cy_beta_alpha2,s,b);

aeroModel.Cz_0_0  = addAdditiveAndMultiplUnc(P_AP2.Cz_0_0,s, b);
aeroModel.Cz_0_alpha  = addAdditiveAndMultiplUnc(P_AP2.Cz_0_alpha,s, b);
aeroModel.Cz_0_alpha2  = addAdditiveAndMultiplUnc(P_AP2.Cz_0_alpha2,s,b);

%% Tether parameter
% Uwes paper
T.c0 =  614600 ;
T.d0 = 473;
T.rho_t = 0.518/170;
T.d_tether = 0.004;
T.CD_tether = 1.2; %AP2 (G.L thesis)
T.np = 5;
T.CD_tether = 1.2; %AP2 (G.L thesis)
T.E = 5.3e9; % from Fagiano
T.eps = 0.02;
T.springconstant_un = T.E*pi*T.d_tether^2/(4*0.02);
T.d_tether = 2e-3; %also AP2 (G.L thesis)

T.l0 = simInit.tether_inital_lenght/(T.np+1);
simInit.pos_p_init = [];
e_t = simInit.pos_W_init/norm(simInit.pos_W_init);
for p = 1 : T.np
    simInit.pos_p_init = [simInit.pos_p_init; p*e_t*T.l0];
end
simInit.vel_p_init = 0 * simInit.pos_p_init;

%% Winch parameters
winchParameter.radius = 0.1;
winchParameter.inertia = 0.08;
winchParameter.friction = 0.6;
simInit.winch_angle_init = simInit.tether_inital_lenght/winchParameter.radius;

%% ================ Control parameters ================
params.mu_a_cmd_filter_traction = 10; 
params.mu_a_cmd_filter_retraction = 10; 
params.alpha_a_cmd_filter_traction = 10; 
params.alpha_a_cmd_filter_retraction = 10; 

params.w0_mu_a_emulation = 3; 
params.w0_alpha_a_emulation = 3; 

params.w0_chi_retraction_max= 1.5000;
params.w0_gamma_retraction_max= 1.0041;
params.w0_chi_retraction_min= 0.3014;
params.w0_gamma_retraction_min= 0.4016;
params.Kp_chi_tau_traction= 0.1867;
params.Kp_gamma_tau_traction= 0.0450;
params.Ki_gamma_tau_traction= 0.0043;
params.Kp_chi_retraction= 2.4000;
params.Kp_gamma_retraction= 2;
params.Ki_chi_retraction= 0.0100;
params.Ki_gamma_retraction= 0.1000;
params.kp_winch= 0.9660;
params.ki_winch= 0.0258;

params.s_retrac_trigger = pi/2;

params.l_tether_max= 700;
params.l_tether_min= 300;

% Guidance traction phase
params.kp_delta_course= 0.0500;
params.kp_delta_gamma= 0.0600;

% Guidance traction phase
params.delta0_offset= 0.15;%0.1500;
params.delta0_slope= 0;

params.Kp_chi_k_traction= 1;
params.Ki_chi_k_traction= 0.1000;

params.Kp_gamma_k_traction= 1;
params.Ki_gamma_k_traction= 0.1000;

% Retraction phase 
params.vamin_max = 30; 
params.vamin = 25; 
params.vamax = 45; 
params.vamax_min = 30; 
params.KpVaControl = 3; 

% maximum shortest tether length
params.l_t_min_max = 300; 

%% Initializing attitudes etc 
% Initializing
M_OW = [cos(ENVMT.windDirection_rad), sin(ENVMT.windDirection_rad), 0;
    sin(ENVMT.windDirection_rad), -cos(ENVMT.windDirection_rad), 0;
    0, 0, -1];
M_tauW = [-sin(lat_in)*cos(x0_sim(10)), -sin(lat_in)*sin(x0_sim(10)), cos(lat_in);
    -sin(x0_sim(10)), cos(x0_sim(10)), 0;
    -cos(lat_in)*cos(x0_sim(10)), -cos(lat_in)*sin(x0_sim(10)), -sin(lat_in)];
phi_t = x0_sim(4); theta_t = x0_sim(5); psi_t = x0_sim(6);
M_tauB = [ cos(psi_t)*cos(theta_t), cos(psi_t)*sin(phi_t)*sin(theta_t) - cos(phi_t)*sin(psi_t), sin(phi_t)*sin(psi_t) + cos(phi_t)*cos(psi_t)*sin(theta_t);
    cos(theta_t)*sin(psi_t), cos(phi_t)*cos(psi_t) + sin(phi_t)*sin(psi_t)*sin(theta_t), cos(phi_t)*sin(psi_t)*sin(theta_t) - cos(psi_t)*sin(phi_t);
    -sin(theta_t),                              cos(theta_t)*sin(phi_t),                              cos(phi_t)*cos(theta_t)];
simInit.M_OB = M_OW * M_tauW' * M_tauB;
M_BO = simInit.M_OB';
phi = atan2( simInit.M_OB(3,2) , simInit.M_OB(3,3) );
theta = -asin( simInit.M_OB(3,1) );
psi = atan2( simInit.M_OB(2,1) , simInit.M_OB(1,1) );
simInit.EULER_init = [phi;theta;psi];
vel_w_W = [12;0;0];
vel_w_B = M_BO*M_OW*vel_w_W;
alpha = x0_sim(3);
beta = x0_sim(2);
Va = x0_sim(1);
M_AB = [cos(alpha)*cos(beta), sin(beta), sin(alpha)*cos(beta);
    -cos(alpha)*sin(beta), cos(beta), -sin(alpha)*sin(beta);
    -sin(alpha), 0, cos(alpha)];
M_BA = M_AB';
simInit.v_k_B_init = vel_w_B + M_BA*[Va;0;0];
simInit.v_k_O_init = simInit.M_OB*simInit.v_k_B_init;
v_a_O = simInit.v_k_O_init - transformFromWtoO( ENVMT.windDirection_rad, vel_w_W );
simInit.v_k_W_init = M_OW'*simInit.v_k_O_init;
simInit.pathangles_k_init(1) = atan2( simInit.v_k_O_init(2), simInit.v_k_O_init(1) );
simInit.pathangles_k_init(2) = -asin( simInit.v_k_O_init(3)/norm(simInit.v_k_O_init) );
chi_a = atan2( v_a_O(2), v_a_O(1) );
gamma_a = -asin( v_a_O(3)/norm(v_a_O) );
simInit.pathangles_a = [chi_a; gamma_a];
simInit.mu_a_init = asin( max( min( (cos(simInit.EULER_init(2))*sin(simInit.EULER_init(1)))/(cos(gamma_a)*cos(0)) + tan(gamma_a) * tan(0), 1),-1) );


end