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
% :Author: Sebastian Rapp (s.rapp@tudelft.nl)

function [pos_O] = transformFromWtoO(windDirection_rad,pos_W)
%TRANSFORMFROMOTOW Summary of this function goes here
%   Detailed explanation goes here

M_OW = [cos(windDirection_rad), sin(windDirection_rad), 0;
        sin(windDirection_rad), -cos(windDirection_rad), 0;
        0, 0, -1];

pos_O = M_OW*pos_W;

end
