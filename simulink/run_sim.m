% Copyright (C) 2023  Nikolaus Vertovec
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
% :Revision: 31-January-2021
% :Author: Nikolaus Vertovec (nikolaus.vertovec@eng.ox.ac.uk)
% :Adapted from: Sebastian Rapp (s.rapp@tudelft.nl) and Dylan Eijkelhof (d.eijkelhof@tudelft.nl)

rmpath('../hj_reachability')
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

Ft_max = 1870;
params.v_reel_retract = constr.winchParameter.v_min;
safety_controller = false;
stopOnRupture = false;

turbulence = true;

% Load data
currentFile = mfilename( 'fullpath' );
[pathstr,~,~] = fileparts( currentFile );
 
load([pathstr, '/../hj_reachability/awe/tables.mat'])
 
random_seed = [ 24335       17059       17670        1747]; %random_seed_mat(idx,:); 
save('random_seed.mat', 'random_seed'); 
%% Define Drag coefficients based on Rapp 2019
P_AP2.initAircraft.alpha          = [-5     0       5       10      15      20      25      30      35];
P_AP2.initAircraft.wing_cL_Static = [0.3    0.8     1.2     1.5     1.5     1.5     1.5     1.5     1.4];
P_AP2.initAircraft.wing_cD_Static = [0.05   0.6     0.1     0.15    2.7     0.42    0.65    0.85    1];
%% Ready to go...
idx =3;  
simInit.perturbed_aero_flag = 0;
simInit.perturbed_ft_moment_fl0ag = 0; 

simOut = sim(name_simulink_model,'ReturnWorkspaceOutputs','on', 'CaptureErrors', 'On');

save(['yout_test_',name_simulink_model,'.mat'], 'simOut');
if make_video
    addpath(genpath('video')); 
    makevideo
end
%% Plotting
number_of_clycles = 3;
if ~isinf(number_of_clycles) % && ~simOut.rupture.Data
    % Power and flight path for last pumping cycle
    P_mech_last_cycle = extractSignalOfLastCycle2(simOut.P_mech, ...
        unique(simOut.cycle_signal_counter.Data), simInit, number_of_clycles);
    Path_last_cycle = extractSignalOfLastCycle_nD(simOut.pos_O, ...
        unique(simOut.cycle_signal_counter.Data), simInit, number_of_clycles);

    Average_power = mean(P_mech_last_cycle.Data);
    [fig_PO, fig_flightpower] = Offline_visualisation_power(...
        P_mech_last_cycle,Path_last_cycle);
    fig_PO = Theoretical_Pcheck(simOut,constr,ENVMT,P_AP2,simInit,T,params,fig_PO, number_of_clycles);
    
    if safety_controller
        % Safety switching for last pumping cycle
        safety_last_cycle = extractSignalOfLastCycle2(simOut.safety_on, ...
            unique(simOut.cycle_signal_counter.Data), simInit, number_of_clycles);
        [fig_saf] = Offline_visualisation_safety(...
            safety_last_cycle,Path_last_cycle);
    end

    % Force
    force_last_cycle = extractSignalOfLastCycle2(simOut.TetherForce, ...
        unique(simOut.cycle_signal_counter.Data), simInit, number_of_clycles);
    [figF, fig_force] = Offline_visualisation_force(...
        force_last_cycle,Path_last_cycle, Ft_max);
else
    if ~simOut.rupture.Data
        % Power and flight path
        Average_power = mean(simOut.P_mech.Data);
        [fig_PO, fig_flightpower] = Offline_visualisation_power(...
            simOut.P_mech,simOut.pos_O);
        fig_PO = Theoretical_Pcheck(simOut,constr,ENVMT,P_AP2,simInit,T,params,fig_PO);
    end

    if safety_controller
        % Safety switching
        [fig_saf] = Offline_visualisation_safety(...
            simOut.safety_on,simOut.pos_O);
    end

    % Force
    [figF, fig_force] = Offline_visualisation_force(...
        simOut.TetherForce,simOut.pos_O, Ft_max, simOut.rupture.Data);
end