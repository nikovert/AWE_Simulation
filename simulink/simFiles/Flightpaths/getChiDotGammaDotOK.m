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
% :Author: Sebastian Rapp (s.rapp@tudelft.nl)

function [chi_dot, gamma_dot] = getChiDotGammaDotOK(windDirection_rad,lat,long, v_k_W, chi_tau_dot, r, chi_k, gamma_k, v_ideal, pos_W, gamma_tau,gamma_tau_dot)

v_k_W  = v_k_W / norm(pos_W); 

v_k_O = transformFromWtoO(windDirection_rad, v_k_W); 
ex = v_k_O/norm(v_k_O); 
pos_O = transformFromWtoO(windDirection_rad, pos_W); 
ez_tmp = -pos_O/norm(pos_O);
ey = cross(ez_tmp, ex); ey = ey/norm(ey);
ez = cross(ex, ey);

M_KbarO = [ex';ey';ez'];    
M_OKbar = M_KbarO'; 

omega_OKbar_Kbar = [-chi_tau_dot*sin(gamma_tau);
                     gamma_tau_dot;
                     chi_tau_dot*cos(gamma_tau)];

v_k_tau = transformFromWtoTau( long, lat, v_ideal); 

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
end



