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

function [x, y, z] = sph2cart(long, lat, r)
%SPH2CART Transform spherical to Cartesian coordinates.
%   [X,Y,Z] = SPH2CART(TH,PSI,R) transforms corresponding elements of
%   data stored in spherical coordinates (azimuth TH, inclination PSI,
%   radius R) to Cartesian coordinates X,Y,Z.  The arrays TH, PSI, and
%   R must be the same size (or any of them can be scalar).  TH and
%   PHI must be in radians.
%
%   TH is the counterclockwise angle in the xy plane measured from the
%   positive x axis.  PSI is the inclination angle from the positive x
%   axis.
%
%   Class support for inputs TH,PHI,R:
%      float: double, single
%
%   See also CART2SPH

x = r .* cos(lat) .* cos(long);
y = r .* cos(lat) .* sin(long);
z = r .* sin(lat);
end