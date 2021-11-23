function [ d_segment ] = calcSegmentDrag( p_d, v_a_p, rho_air, CD_tether,d_tether )
%CALCSEGMENTDRAG Summary of this function goes here
%   Detailed explanation goes here

%% Readme
% Equations are adapted from Fechner et al, Dynamic Model of a Pumping Kite
% Power System, Renewable Energy, 2015.
% Implementation: Sebastian Rapp, Wind Energy Institute, Faculty of
% Aerospace Engineering, TU Delft
% Mail: s.rapp@tudelft.nl
% Last change: 19.01.2018
%
% General description:
% This functions calculates the segment drag.
% 
% Inputs:
% CD_tether: Tether drag coefficient
% d_tether: Tether diameter
% rho_air: air density
% p_d: relative postion vector between two particles
% v_a_p: relative airspeed at one particle
%
% Outputs:
% d_segment: segment drag


dp_norm = p_d - p_d .* v_a_p .* v_a_p / (v_a_p' * v_a_p) ; % projection of segment direction perpedicular to v_a direction
va_perp = v_a_p - ( p_d' *v_a_p )/norm(p_d) * p_d/norm(p_d); 
A_eff = norm(dp_norm) * d_tether; % projected tether area perpendicular to v_a
d_segment = 0.5 * rho_air * CD_tether * va_perp * norm(va_perp) * A_eff; % particle drag


end


