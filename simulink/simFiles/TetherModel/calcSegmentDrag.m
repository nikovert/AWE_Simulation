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

function [ d_segment ] = calcSegmentDrag( p_d, v_a_p, rho_air, CD_tether,d_tether )
%CALCSEGMENTDRAG Summary of this function goes here
%   Detailed explanation goes here

%% Readme
% Equations are adapted from Fechner et al, Dynamic Model of a Pumping Kite
% Power System, Renewable Energy, 2015.
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


