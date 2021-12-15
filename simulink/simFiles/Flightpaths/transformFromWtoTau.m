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

function vec_tau = transformFromWtoTau(  lamb,phi, vec_W )

M_tauW = [-sin(phi)*cos(lamb), -sin(phi)*sin(lamb), cos(phi);
    -sin(lamb), cos(lamb), 0;
    -cos(phi)*cos(lamb), -cos(phi)*sin(lamb), -sin(phi)];

vec_tau = M_tauW * vec_W; 