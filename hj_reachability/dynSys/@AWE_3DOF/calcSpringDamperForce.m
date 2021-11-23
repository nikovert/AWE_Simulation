function [ fsd_tot_vec_i ] = calcSpringDamperForce(obj, c, p_norm, tether_diff, d, p_d, v_d )
%CALC_SPRINGDAMPERFORCE calcSpringDamperForce( c, p_norm, l_s, d, p_d, v_d )

% Equations are adapted from Fechner et al, Dynamic Model of a Pumping Kite
% Power System, Renewable Energy, 2015.
% Implementation: Sebastian Rapp, Wind Energy Institute, Faculty of
% Aerospace Engineering, TU Delft
% Mail: s.rapp@tudelft.nl
% Last change: 19.01.2018
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
