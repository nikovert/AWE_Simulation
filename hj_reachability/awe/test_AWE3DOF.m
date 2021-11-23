%% Add Path
addpath('../dynSys')
%% Initialize
Va      = 31.3466;
alpha   = 0.1222;
mu_a    = 0.0834;
chi_a   = -2.2106;
gamma_a = -0.1147;
long    = -9.5681e-17;
lat     = 1.3963;
h_tau   = 250;

x = [long, lat, h_tau, Va, chi_a, gamma_a]';
u = [alpha, mu_a]';

dynSys = AWE_3DOF(x);
dynSys.v_w_O = [-15.8031   -0.0001   -0.0000]';
x_dot = dynSys.dynamics(0, x, u, 0);
