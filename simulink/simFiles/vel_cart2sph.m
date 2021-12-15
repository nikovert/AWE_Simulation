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
% :Author: Nikolaus Vertovec (nikolaus.vertovec@eng.ox.ac.uk)

function [long_dot, lat_dot, r_dot] = vel_cart2sph(x, y, z, vx, vy, vz)
% VEL_CART2SPH transforms cartesian velocity to spherical coordinates
    r_dot    = (x.*vx + y.*vy + z.*vz)./(sqrt(x.^2 + y.^2 + z.^2));
    long_dot = (x.*vy - y.*vx)./(x.^2 + y.^2);
    lat_dot  = -(z.*x.*vx + z.*y.*vy - (x.^2 + y.^2).*vz)./((x.^2 + y.^2 + z.^2).^(3./2) .* sqrt(1 - z.^2 ./ (x.^2 + y.^2 + z.^2)));
end