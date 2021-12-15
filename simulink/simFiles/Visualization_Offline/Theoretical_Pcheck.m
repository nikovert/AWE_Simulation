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
% :Adapted from: Dylan Eijkelhof (d.eijkelhof@tudelft.nl)

function fig_PO = Theoretical_Pcheck(simOut,constr,ENVMT,P_AP2,simInit,T,params,fig_PO ,number_of_clycles)
% Theoretical check, add values to instantaneous power plot.
%
% :param simOut: Simulation output
% :param constr: Aircraft manoeuvre and winch constraints
% :param ENVMT: Environmental parameters
% :param P_AP2: Aircraft parameters
% :param simInit: Simulation initialisation parameters
% :param T: Tether dimensions and material properties
% :param params: Flight/Winch controller parameters
% :param fig_PO: Instantaneous power figure (Offline_visualisation_power.m)
% :param number_of_clycles: the umber of sampels cycles considered
% :returns: 
%           - **fig_PO** - Instantaneous power plot of converged power cycle with theoretical power.
%
% .. note::
%           - Loyd peak power with cosine losses from cycle CL and CD, 
%           - Costello et al. average power from cycle CL and CD

%------------- BEGIN CODE --------------
if nargin<8
    fig_PO = NaN;
end

if nargin < 4
    number_of_clycles = 1;
end

%% Calculate parameters needed to fill in the theoretical equation of Loyd
colourvec = {'#77AC30', '#A2142F', '#EDB120', '#0072BD', 'm', 'c', 'y'};
TL = extractSignalOfLastCycle_nD(simOut.tether_length, ...
            unique(simOut.cycle_signal_counter.Data), simInit, number_of_clycles);
Vw = extractSignalOfLastCycle_nD(simOut.v_w, ...
        unique(simOut.cycle_signal_counter.Data), simInit, number_of_clycles);
Flightstate_last_cycle = extractSignalOfLastCycle_nD(simOut.sub_flight_state, ...
        unique(simOut.cycle_signal_counter.Data), simInit, number_of_clycles);
Alpha_last_cycle = extractSignalOfLastCycle_nD(simOut.alpha, ...
        unique(simOut.cycle_signal_counter.Data), simInit, number_of_clycles);
% Lift_last_cycle = extractSignalOfLastCycle_nD(simOut.lift, ...
%         unique(simOut.cycle_signal_counter.Data), simInit, number_of_clycles);
Drag_last_cycle = extractSignalOfLastCycle_nD(simOut.drag, ...
        unique(simOut.cycle_signal_counter.Data), simInit, number_of_clycles);
Faero_last_cycle = extractSignalOfLastCycle_nD(simOut.Faero, ...
        unique(simOut.cycle_signal_counter.Data), simInit, number_of_clycles);
P_mech_last_cycle = extractSignalOfLastCycle_nD(simOut.P_mech, ...
        unique(simOut.cycle_signal_counter.Data), simInit, number_of_clycles);
Phi_last_cycle = extractSignalOfLastCycle_nD(simOut.phi0_mean_path, ...
        unique(simOut.cycle_signal_counter.Data), simInit, number_of_clycles);
dp_retract = find((Flightstate_last_cycle.Data==6|Flightstate_last_cycle.Data==7),1);
t_retract = Flightstate_last_cycle.Time(dp_retract);    

amin = rad2deg(constr.alpha_a_min);
amax = rad2deg(constr.alpha_a_max);
ind = (P_AP2.initAircraft.alpha>=amin & P_AP2.initAircraft.alpha<=amax);
x = P_AP2.initAircraft.alpha(ind);
y1 = P_AP2.initAircraft.wing_cL_Static(ind);
y2 = P_AP2.initAircraft.wing_cD_Static(ind);
alpha_list = rad2deg(Alpha_last_cycle.Data(Alpha_last_cycle.Time<=t_retract));
alpha_list(alpha_list>amax) = amax;
alpha_list(alpha_list<amin) = amin;
C_L = interp1(x,y1,alpha_list);
C_D = interp1(x,y2,alpha_list);
C_D_tether = T.CD_tether*TL.Data(TL.Time<=t_retract)*T.d_tether/(4*P_AP2.S_wing);
         
[cl3cd2,itether] = max(C_L.^3./(C_D+C_D_tether).^2);

% Just reel-out phase is taken into account as this is where peak power is
% expected to occur.
%p_loyd_tether = 2/27*ENVMT.rhos*P_AP2.S_wing*mean(Vw.Data)^3*cl3cd2*cos(params.phi0_booth)^3/1e3;

disp(['CL @ CDeff: ', num2str(C_L(itether))])
disp(['CDeff: ', num2str(C_D(itether)+C_D_tether(itether))])
tltest = TL.Data(TL.Time<=t_retract);
disp(['Tether length: ', num2str(tltest(itether))])
%disp(['Peak power (reel-out): ', num2str(p_loyd_tether), ' kW'])
Pw = 0.5*ENVMT.rhos*mean(Vw.Data)^3;
disp(['Power density: ', num2str(Pw)])

%% Costello et al.
% Lift = Lift_last_cycle.Data(Lift_last_cycle.Time<=t_retract);
Drag = (Drag_last_cycle.Data(Drag_last_cycle.Time<=t_retract));
Faero = (Faero_last_cycle.Data(Faero_last_cycle.Time<=t_retract));
v_wind = (Vw.Data(Vw.Time<=t_retract));
Phi = (Phi_last_cycle.Data(Phi_last_cycle.Time<=t_retract));
Pw = 0.5*ENVMT.rhos*(v_wind).^3;
e = (cos(Phi+asin(...
    (Drag./Faero).*sin(Phi)+...
    (P_AP2.mass*ENVMT.g./Faero).*cos(Phi))).^3);
zeta = 4/27*(C_L.^3./(C_D+C_D_tether).^2);
Costello_Pavg_max = mean(e.*Pw.*P_AP2.S_wing.*zeta/1e3);
disp(['Costello Pavg max. (reel-out): ', num2str(Costello_Pavg_max), ' kW'])
% Pout = mean([P_mech_last_cycle.Data(1:dp_retract-1);0.*P_mech_last_cycle.Data(dp_retract:end)]);

disp(['zeta: ', num2str(mean(zeta))])
disp(['e: ', num2str(mean(e))])

%% Plotting
if ishandle(fig_PO)
    figure(fig_PO)
    legend('Interpreter', 'latex')
    
    if Costello_Pavg_max < 1
        mpc = [num2str(Costello_Pavg_max*1e3,3) ' W'];
    else
        mpc = [num2str(Costello_Pavg_max,3) ' kW'];
    end
    plot([0, Vw.Time(end)-Vw.Time(1)],[Costello_Pavg_max,Costello_Pavg_max],'-.','Color',colourvec{5},'DisplayName',['$\tilde{P}$ (' mpc ')'])
    
end
end