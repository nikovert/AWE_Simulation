function [long_dot, lat_dot, r_dot] = vel_cart2sph(x, y, z, vx, vy, vz)
    r_dot    = (x.*vx + y.*vy + z.*vz)./(sqrt(x.^2 + y.^2 + z.^2));
    long_dot = (x.*vy - y.*vx)./(x.^2 + y.^2);
    lat_dot  = -(z.*x.*vx + z.*y.*vy - (x.^2 + y.^2).*vz)./((x.^2 + y.^2 + z.^2).^(3./2) .* sqrt(1 - z.^2 ./ (x.^2 + y.^2 + z.^2)));
end