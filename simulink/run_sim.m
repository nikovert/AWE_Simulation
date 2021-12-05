clear 
close all
clc
%%
% -------------------------------------------Setup -------------------------------------------
pointmass_sim_flag = 1;
complex_tether_flag = 1; % only for the PM sim  
controller_flag = 1;
make_video = false;

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

turbulence = true;
if turbulence
    params.Ft_set_traction = 1500;
end

% Load data
currentFile = mfilename( 'fullpath' );
[pathstr,~,~] = fileparts( currentFile );
 
load([pathstr, '/../hj_reachability/awe/tables.mat'])
 
random_seed = [ 24335       17059       17670        1747]; %random_seed_mat(idx,:); 
save('random_seed.mat', 'random_seed'); 
%% Ready to go...
idx =3;  
simInit.perturbed_aero_flag = 0;
simInit.perturbed_ft_moment_flag = 0; 

simOut = sim(name_simulink_model,'ReturnWorkspaceOutputs','on', 'CaptureErrors', 'On');

save(['yout_test_',name_simulink_model,'.mat'], 'simOut');
if make_video
    addpath(genpath('video')); 
    makevideo
end

%% Power and flight path for last pumping cycle
P_mech_last_cycle = extractSignalOfLastCycle2(simOut.P_mech, ...
    simOut.cycle_signal_counter, simInit);
Average_power = mean(P_mech_last_cycle.Data);

Path_last_cycle = extractSignalOfLastCycle_nD(simOut.pos_O, ...
    simOut.cycle_signal_counter, simInit );
[fig_PO, fig_flightpower] = Offline_visualisation_power(...
            P_mech_last_cycle,Path_last_cycle);
%fig_PO = Theoretical_Pcheck(simOut,constr,ENVMT,P_AP2,simInit,T,params,fig_PO);