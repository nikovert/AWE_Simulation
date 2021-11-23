function [x, y, z] = sph2cart(long, lat, r)
%SPH2CART Transform spherical to Cartesian coordinates.
%   [X,Y,Z] = SPH2CART(TH,PSI,R) transforms corresponding elements of
%   data stored in spherical coordinates (azimuth TH, inclination PSI,
%   radius R) to Cartesian coordinates X,Y,Z.  The arrays TH, PSI, and
%   R must be the same size (or any of them can be scalar).  TH and
%   PHI must be in radians.
%
%   TH is the counterclockwise angle in the xy plane measured from the
%   positive x axis.  PSI is the inclination angle from the positive z
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