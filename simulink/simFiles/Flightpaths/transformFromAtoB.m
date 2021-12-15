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

function vec_B = transformFromAtoB( alpha, beta,  vec_A )

c_a = cos(alpha);
s_a = sin(alpha);

c_b = cos(beta);
s_b = sin(beta);

M_A2B = zeros(3,3);

M_A2B(1,1) = c_a*c_b;
M_A2B(1,2) = -c_a*s_b;
M_A2B(1,3) = -s_a;
M_A2B(2,1) = s_b;
M_A2B(2,2) = c_b;
M_A2B(2,3) = 0;
M_A2B(3,1) = s_a*c_b;
M_A2B(3,2) = -s_a*s_b;
M_A2B(3,3) = c_a;

vec_B = M_A2B * vec_A; 