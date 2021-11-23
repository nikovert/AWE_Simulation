function [pos_O] = transformFromWtoO(windDirection_rad,pos_W)
%TRANSFORMFROMOTOW Summary of this function goes here
%   Detailed explanation goes here

M_OW = [cos(windDirection_rad), sin(windDirection_rad), 0;
        sin(windDirection_rad), -cos(windDirection_rad), 0;
        0, 0, -1];

pos_O = M_OW*pos_W;

end
