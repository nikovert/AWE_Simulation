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

function [vx, vy, vz] = vel_sph2cart(long, lat, r, long_dot,lat_dot,r_dot)
%SPH2CART Transform spherical velocities to Cartesian coordinates.

vx = r_dot .* cos(lat) .* cos(long) - r .* sin(lat) .* lat_dot .* cos(long) - r .* cos(lat) .* sin(long) .* long_dot;
vy = r_dot .* cos(lat) .* sin(long) - r .* sin(lat) .* lat_dot .* sin(long) + r .* cos(lat) .* cos(long) .* long_dot;
vz = r_dot .* sin(lat) + r .* cos(lat) .* lat_dot;
end