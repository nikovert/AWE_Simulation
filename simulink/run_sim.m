clear 
close all
clc
%%
% -------------------------------------------Setup -------------------------------------------
pointmass_sim_flag = 1;
complex_tether_flag = 1; % only for the PM sim  
controller_flag = 1;
delete('log_*'); 
delete('yout*'); 
addpath(genpath('simFiles')); 

name_simulink_model = 'AWES3_CL_wTether_v4_PM';
lat_init = 80*pi/180; 

% Initialize parameters for simulation
[x0_sim,u0_sim, act, aeroModel, base_windspeed, constr,...
    ENVMT, P_AP2, simInit, T, winchParameter,params] = initAllSimParams(lat_init);
params.direction = -1;

Ft_max = 1865;
safety_controller = true;

% Load data
currentFile = mfilename( 'fullpath' );
[pathstr,~,~] = fileparts( currentFile );
 
load([pathstr, '/../HJReachability/awe/tables.mat'])
 
%% Ready to go...
idx =3;  
simInit.perturbed_aero_flag = 0;
simInit.perturbed_ft_moment_flag = 0; 

yout = sim(name_simulink_model,'ReturnWorkspaceOutputs','on', 'CaptureErrors', 'On');

save(['yout_test_',name_simulink_model,'.mat'], 'yout');

addpath(genpath('video')); 
makevideo

delete(gcp('nocreate'));
