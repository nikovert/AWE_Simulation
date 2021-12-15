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

function e_course = wrapCourseError(chi_ref, chi)
%#codegen
% The vectors are lying in the tangential plane attached to the current
% position on the small earth.
p_bearing_vec =   [cos(chi);         sin(chi);         0];
ref_bearing_vec = [cos(chi_ref); sin(chi_ref); 0];

e_z = cross( ref_bearing_vec, p_bearing_vec ); 

dot_product = dot( ref_bearing_vec, p_bearing_vec ); 
if dot_product > 1 
    dot_product = 1; 
end
if dot_product < -1 
    dot_product = -1; 
end
e_course = - sign( e_z(3) ) * acos( dot_product );
