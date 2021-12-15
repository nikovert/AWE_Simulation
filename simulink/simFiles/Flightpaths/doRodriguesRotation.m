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

function vrot = doRodriguesRotation(pos_W, p_target_W, v)

theta = acos( max( min( pos_W'*p_target_W / norm(pos_W) / norm( p_target_W ), 1),-1) );
k = cross( p_target_W, pos_W );

if abs(theta) < 1e-12 || norm(k) < 1e-12
    vrot = v;
else
    k = k/norm(k);
    vrot = v*cos(theta)+cross(k,v)*sin(theta)+k * (k'*v)*(1-cos(theta));
end
end