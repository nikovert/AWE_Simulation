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
% :Adapted from: Sebastian Rapp (s.rapp@tudelft.nl)

function [ fsd_tot_vec_i ] = calcSpringDamperForce( c, p_norm, l_s, d, p_d, v_d )
%CALC_SPRINGDAMPERFORCE Summary of this function goes here

%% Readme
% Equations are adapted from Fechner et al, Dynamic Model of a Pumping Kite
% Power System, Renewable Energy, 2015.
%
% General description:
% This functions calculates the spring-damper force acting in one tether segment.
% 
% Inputs:
% c: spring stiffness
% p_norm: euclidean distance between two adjacent tether particles
% l_s: segment length
% d: damping coefficient
% p_d: relative postion vector between two particles
% v_d: L2-norm velocities of two adjacent particles
%
% Outputs:
% Spring-Damper-force

x  =  p_norm - l_s;
epsilon = 0; 
delta = 1e-3; 

if x<=0 % no compressive forces allowed
    c = c * min( max( 1/delta * ( x + epsilon )+1, 0), 1); % smooth
    d = d * min( max( 1/delta * ( x + epsilon )+1, 0), 1); 
    fsd_tot_vec_i = (c * ( x ) + d * ( p_d/p_norm)' * v_d ) * p_d/p_norm;
else
    fsd_tot_vec_i = (c * ( x ) + d * ( p_d/p_norm)' * v_d ) * p_d/p_norm;
end

end
