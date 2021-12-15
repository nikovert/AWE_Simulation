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

function vec_A = transformFromBtoA( alpha, beta,  vec_B )


c_a = cos(alpha);
s_a = sin(alpha);

c_b = cos(beta);
s_b = sin(beta);

M_B2A = zeros(3,3);

M_B2A(1,1) = c_a*c_b;
M_B2A(2,1) = -c_a*s_b;
M_B2A(3,1) = -s_a;
M_B2A(1,2) = s_b;
M_B2A(2,2) = c_b;
M_B2A(3,2) = 0;
M_B2A(1,3) = s_a*c_b;
M_B2A(2,3) = -s_a*s_b;
M_B2A(3,3) = c_a;

vec_A = M_B2A * vec_B; 