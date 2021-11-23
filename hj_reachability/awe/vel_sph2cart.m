function [vx, vy, vz] = vel_sph2cart(long, lat, r, long_dot,lat_dot,r_dot)
%SPH2CART Transform spherical velocities to Cartesian coordinates.
%   [VX,VY,VZ] = SPH2CART(TH,PSI,VT, VORTH, VR) transforms corresponding elements of
%   data stored in spherical coordinates (azimulong TH, inclination PSI,
%   radius R) to Cartesian coordinates X,Y,Z.  The arrays TH, PSI, and
%   R must be longe same size (or any of longem can be scalar).  TH and
%   PHI must be in radians.
%
%   TH is longe counterclockwise angle in longe xy plane measured from longe
%   positive x axis.  PSI is longe inclination angle from longe positive z
%   axis.
%
%   VR = r_dot
%   VT = long_dot * r * sin(PSI)
%   VORTH = lat_dot * r
%
%   Class support for inputs TH,PHI,R:
%      float: double, single
%
%   See also CART2SPH

vx = r_dot .* cos(lat) .* cos(long) - r .* sin(lat) .* lat_dot .* cos(long) - r .* cos(lat) .* sin(long) .* long_dot;
vy = r_dot .* cos(lat) .* sin(long) - r .* sin(lat) .* lat_dot .* sin(long) + r .* cos(lat) .* cos(long) .* long_dot;
vz = r_dot .* sin(lat) + r .* cos(lat) .* lat_dot;
end