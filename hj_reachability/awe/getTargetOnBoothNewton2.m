function [ sol,p_C_P ,nonconverged_points_list] = getTargetOnBoothNewton2(Lem, pos_P, s_old, distanceOrigin, direction)
    [sol,exceedMaxIter, nonconverged_points_list] = doNewtonIterationBooth2(s_old,Lem, pos_P, direction);

    long_path = Lem.b .* sin(sol) ./ ( 1+(Lem.a./Lem.b .* cos(sol)).^2 );
    lat_path  = Lem.a .* sin(sol) .* cos(sol) ./ ( 1+(Lem.a./Lem.b .* cos(sol)).^2 ) ;
    
    [pos_C_P_x,pos_C_P_y,pos_C_P_z] = sph2cart(long_path,lat_path,distanceOrigin);
    p_C_P = {pos_C_P_x;pos_C_P_y;pos_C_P_z};
end