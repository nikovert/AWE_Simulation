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

function [Chi_cmd,Delta_chi,theta, Delta_chi_dot, deltaDot1] = calculateCommandedCourse2(t,delta_vec, delta, delta0, chi_parallel, L_sol, pos_W, v_w_vec, chi_tau_is)
%CALCULATEDESIREDCOURSE Summary of this function goes here
%   Detailed explanation goes here

% check sign with curve frame - y component is negative in this frame -
% delta Chi is positive and vice versa.
% add a predictive part to smooth highly curved paths.
% basis:
e1 = t/norm(t);
e3 = -L_sol/norm(L_sol);
e2 = cross(e3, e1);
M_CW = [e1';e2';e3'];
pos_C = M_CW * (pos_W-L_sol);

% Law 1: 
Delta_chi = atan2( -sign(pos_C(2)) * delta, delta0 );

% Law 2: 
%Delta_chi = asin( -sign(pos_C(2))*max( min( delta*delta0, 1), -1) ); 
%1;


Chi_cmd = chi_parallel + Delta_chi;

if Chi_cmd > pi
    Chi_cmd = -pi + mod( Chi_cmd, pi );
elseif Chi_cmd < -pi
    Chi_cmd = pi + mod(Chi_cmd, -pi);
end

if sign(pos_C(2)) < 0
    theta = pi/2 - Delta_chi;
else
    theta = pi/2 + Delta_chi;
end
%if delta < 1e-6 
  %  Delta_chi_dot = 0; 
   % deltaDot1 = 0; 
%else  
    %Delta_chi_dot = sign(pos_C(2)) * norm(v_w_vec/norm(pos_W)) / delta0^2 * delta/(1+(delta/delta0)^2)^(3/2);
    %deltaDot1 = sign(pos_C(2)) * norm(v_w_vec/norm(pos_W)) * sin(Delta_chi);
    %deltaDot1 =  (-v_w_vec'/norm(pos_W) * delta_vec );
    e_chi = wrapCourseError( Chi_cmd, chi_tau_is);
    deltaDot1 = -norm( v_w_vec'/norm(pos_W) ) * sign(pos_C(2)) * sin( -Delta_chi + e_chi )  ;
    Delta_chi_dot = -sign(pos_C(2)) / delta0 / (1+(delta/delta0)^2) * deltaDot1; 
%end

end

