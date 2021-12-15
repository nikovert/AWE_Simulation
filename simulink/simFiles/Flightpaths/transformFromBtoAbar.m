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

function vec_Abar = transformFromBtoAbar( mu_a, alpha, beta,  vec_B )

c_m = cos(mu_a);
s_m = sin(mu_a);

c_a = cos(alpha);
s_a = sin(alpha);

c_b = cos(beta);
s_b = sin(beta);

M_B2ABar = zeros(3,3);

M_B2ABar(1,1) = c_a*c_b;
M_B2ABar(1,2) = s_b;
M_B2ABar(1,3) = s_a*c_b;
M_B2ABar(2,1) = -c_a*s_b*c_m+s_a*s_m;
M_B2ABar(2,2) = c_b*c_m;
M_B2ABar(2,3) = -s_a*s_b*c_m-s_m*c_a;
M_B2ABar(3,1) = -c_a*s_b*s_m-s_a*c_m;
M_B2ABar(3,2) = c_b*s_m;
M_B2ABar(3,3) = -s_a*s_b*s_m+c_a*c_m;

vec_Abar = M_B2ABar * vec_B; 