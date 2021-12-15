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

function vec_Abar = transformFromOtoAbar( chi_a,gamma_a,  vec_O )

c_chia = cos(chi_a);
s_chia = sin(chi_a);

c_gama = cos(gamma_a);
s_gama = sin(gamma_a);

M_O2ABar = zeros(3,3);

M_O2ABar(1,1) = c_chia*c_gama;
M_O2ABar(1,2) = s_chia*c_gama;
M_O2ABar(1,3) = -s_gama;
M_O2ABar(2,1) = -s_chia;
M_O2ABar(2,2) = c_chia;
M_O2ABar(2,3) = 0;
M_O2ABar(3,1) = c_chia*s_gama;
M_O2ABar(3,2) = s_chia*s_gama;
M_O2ABar(3,3) = c_gama;

vec_Abar = M_O2ABar * vec_O; 