function [ fsd_tot_vec_i ] = calcSpringDamperForce( c, p_norm, l_s, d, p_d, v_d )
%CALC_SPRINGDAMPERFORCE Summary of this function goes here

%% Readme
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

x  =  p_norm - l_s;
epsilon = 0; 
delta = 1e-3; 

if x<=0 % no compressive forces allowed
    c = c * min( max( 1/delta * ( x + epsilon )+1, 0), 1); % smooth
    d = d * min( max( 1/delta * ( x + epsilon )+1, 0), 1); 
    fsd_tot_vec_i = (c * ( p_norm - l_s ) + d * ( p_d/p_norm)' * v_d ) * p_d/p_norm;
else
    fsd_tot_vec_i = (c * ( p_norm - l_s ) + d * ( p_d/p_norm)' * v_d ) * p_d/p_norm;
end

end
