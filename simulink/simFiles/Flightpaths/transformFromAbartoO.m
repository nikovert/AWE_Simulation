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

function vec_O = transformFromAbartoO( chi_a,gamma_a,  vec_Abar )

c_chia = cos(chi_a);
s_chia = sin(chi_a);

c_gama = cos(gamma_a);
s_gama = sin(gamma_a);

M_ABar2O = zeros(3,3);

M_ABar2O(1,1) = c_chia*c_gama;
M_ABar2O(2,1) = s_chia*c_gama;
M_ABar2O(3,1) = -s_gama;
M_ABar2O(1,2) = -s_chia;
M_ABar2O(2,2) = c_chia;
M_ABar2O(3,2) = 0;
M_ABar2O(1,3) = c_chia*s_gama;
M_ABar2O(2,3) = s_chia*s_gama;
M_ABar2O(3,3) = c_gama;

vec_O = M_ABar2O * vec_Abar; 