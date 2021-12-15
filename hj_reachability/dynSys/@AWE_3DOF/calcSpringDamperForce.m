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

function [ fsd_tot_vec_i ] = calcSpringDamperForce(obj, c, p_norm, tether_diff, d, p_d, v_d )
%CALC_SPRINGDAMPERFORCE calcSpringDamperForce( c, p_norm, l_s, d, p_d, v_d )

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

x  =  tether_diff;
epsilon = 1e-19; 
delta = 1e-3; 


c = (x<=0) .* c .* min( max( 1/delta * ( x + epsilon )+1, 0), 1) + (x>0) .* c; % smooth
d = (x<=0) .* d .* min( max( 1/delta * ( x + epsilon )+1, 0), 1) + (x>0) .* d; 
p_dir = mult_cellMatrix(p_d,{1./(p_norm+epsilon)});
fsd_tot_vec_i = mult_cellMatrix( ...
                    add_cellMatrix( ...
                        {c.*x}, ...
                        mult_cellMatrix({d}, ...
                            mult_cellMatrix(transpose(p_dir), v_d)) ...
                    ), ...
                    p_dir ...
                );


end
