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

function [ d_segment ] = calcSegmentDrag(obj, p_d, v_a_p, rho_air, CD_tether,d_tether )
%CALCSEGMENTDRAG  [ d_segment ] = calcSegmentDrag( p_d, v_a_p, rho_air, CD_tether,d_tether )
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


dp_norm = sub_cellMatrix(p_d,element_div_cellMatrix(element_mult_cellMatrix(element_mult_cellMatrix(p_d,v_a_p),v_a_p), mult_cellMatrix(transpose(v_a_p), v_a_p))); % projection of segment direction perpedicular to v_a direction
va_perp = sub_cellMatrix(v_a_p, ...
            mult_cellMatrix( ...
                element_div_cellMatrix(mult_cellMatrix(transpose(p_d),v_a_p ),{norm_cellVec(p_d)}), ...
                element_div_cellMatrix(p_d,{norm_cellVec(p_d)}) ...
            ) ...
          ); 
A_eff = norm_cellVec(dp_norm) * d_tether; % projected tether area perpendicular to v_a
d_segment = mult_cellMatrix({0.5 * rho_air * CD_tether * norm_cellVec(va_perp) .* A_eff}, va_perp); % particle drag



end


