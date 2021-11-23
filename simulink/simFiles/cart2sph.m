function [long, lat, r] = cart2sph(x, y, z)
%CART2SPH Transform Cartesian to spherical coordinates.
%   [TH,PSI,R] = CART2SPH(X,Y,Z) transforms corresponding elements of
%   data stored in Cartesian coordinates X,Y,Z to spherical
%   coordinates (azimuth TH, inclination PSI, and radius R).  The arrays
%   X,Y, and Z must be the same size (or any of them can be scalar).
%   TH and PSI are returned in radians.
%
%   TH is the counterclockwise angle in the xy plane measured from the
%   positive x axis.  PSI is the inclination angle from the positive z
%   axis.
%
%   Class support for inputs X,Y,Z:
%      float: double, single
%
%   See also SPH2CART

r = sqrt(x.^2+y.^2+z.^2);
long = atan(y./x);
lat = asin(z./r);
end
