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

function vec_W = transformFromTautoW( lamb,phi,  vec_tau )

M_Wtau = [-sin(phi)*cos(lamb),-sin(lamb),-cos(phi)*cos(lamb);
    -sin(phi)*sin(lamb), cos(lamb), -cos(phi)*sin(lamb);
    cos(phi), 0, -sin(phi)];


vec_W = M_Wtau * vec_tau; 