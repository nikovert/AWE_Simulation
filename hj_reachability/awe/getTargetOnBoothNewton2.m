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
% :Adapted from: Sebastian Rapp (s.rapp@tudelft.nl)

function [ sol,p_C_P ,nonconverged_points_list] = getTargetOnBoothNewton2(Lem, pos_P, s_old, distanceOrigin, direction)
    [sol,exceedMaxIter, nonconverged_points_list] = doNewtonIterationBooth2(s_old,Lem, pos_P, direction);

    long_path = Lem.b .* sin(sol) ./ ( 1+(Lem.a./Lem.b .* cos(sol)).^2 );
    lat_path  = Lem.a .* sin(sol) .* cos(sol) ./ ( 1+(Lem.a./Lem.b .* cos(sol)).^2 ) ;
    
    [pos_C_P_x,pos_C_P_y,pos_C_P_z] = sph2cart(long_path,lat_path,distanceOrigin);
    p_C_P = {pos_C_P_x;pos_C_P_y;pos_C_P_z};
end