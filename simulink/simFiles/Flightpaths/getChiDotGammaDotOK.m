function [chi_dot, gamma_dot] = getChiDotGammaDotOK(windDirection_rad,lat,long, v_k_W, chi_tau_dot, r, chi_k, gamma_k, v_ideal, pos_W, gamma_tau,gamma_tau_dot)
%UNTITLED Summary of this function goes here
v_k_W  = v_k_W / norm(pos_W); 
%v_k_tau = transformFromWtoTau( long, lat, v_k_W);
%v_k_tau(3) = 0;
%v_k_W = transformFromTautoW(long,lat, v_k_tau);

v_k_O = transformFromWtoO(windDirection_rad, v_k_W); 
ex = v_k_O/norm(v_k_O); 
pos_O = transformFromWtoO(windDirection_rad, pos_W); 
ez_tmp = -pos_O/norm(pos_O);
ey = cross(ez_tmp, ex); ey = ey/norm(ey);
ez = cross(ex, ey);

% ex_W = transformFromOtoW(windDirection_rad, ex);
% ey_W = transformFromOtoW(windDirection_rad, ey);
% ez_W = transformFromOtoW(windDirection_rad, ez);

M_KbarO = [ex';ey';ez'];    
M_OKbar = M_KbarO'; 

%omega_OKbar_Kbar = [0;0;chi_tau_dot];
%gamma_tau = 0;%5*pi/180;  % gammaDot = 0
omega_OKbar_Kbar = [-chi_tau_dot*sin(gamma_tau);
                     gamma_tau_dot;
                     chi_tau_dot*cos(gamma_tau)];
%omega_OKbar_Kbar = [0;0;chi_tau_dot];

v_k_tau = transformFromWtoTau( long, lat, v_ideal); 
%v_k_tau_is = transformFromWtoTau( long, lat, v_k_W/norm(v_k_W)*norm(v_ideal)); 

long_dot = v_k_tau(2) /(r*cos(lat)); 
lat_dot = v_k_tau(1) / r; 
  
omega_WT_W = [lat_dot*sin(long); 
             -lat_dot*cos(long); 
              long_dot];

omega_WT_O = transformFromWtoO(windDirection_rad, omega_WT_W ); 

omega_OT_O = omega_WT_O; 

omega_OKbar_O = omega_OT_O + M_OKbar * omega_OKbar_Kbar; 
omega_OKbar_Kbar = M_KbarO * omega_OT_O + omega_OKbar_Kbar;

gamma_dot = omega_OKbar_O(1) * (-sin(chi_k)) + omega_OKbar_O(2) * cos(chi_k); 
mu = atan2( M_OKbar(3,2),M_OKbar(3,3) ); 
chi_dot = (omega_OKbar_Kbar(2) * sin(mu) + omega_OKbar_Kbar(3) * cos(mu) ) / cos(gamma_k); 
% compare with omega_OB_B definition with Theta = gamma, Phi = mu, Psi =
% chi.
end



