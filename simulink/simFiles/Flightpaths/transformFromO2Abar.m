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

function vec_Abar = transformFromO2Abar( chi_a, gamma_a, vec_O )

M_AbarO = [cos(chi_a)*cos(gamma_a), sin(chi_a)*cos(gamma_a), -sin(gamma_a); 
    -sin(chi_a), cos(chi_a), 0; 
    cos(chi_a)*sin(gamma_a), sin(chi_a)*sin(gamma_a), cos(gamma_a)];
vec_Abar = M_AbarO*vec_O; 