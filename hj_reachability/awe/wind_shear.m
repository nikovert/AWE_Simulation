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

function v_w_O = wind_shear(h,windDirection_rad, base_windspeed)
%Wind_Shear This implementation is based on the mathematical representation in the Military Specification MIL-F-8785C
% h should be in meters, as it will be converted to feet
z_0 = 0.15;
vw = base_windspeed * log(3.281 * h/z_0)/log(20/z_0);
v_w_O = cell(3,1);
v_w_O{1} = cos(windDirection_rad) * vw;
v_w_O{2} = sin(windDirection_rad) * vw;
v_w_O{3} = 0;
end