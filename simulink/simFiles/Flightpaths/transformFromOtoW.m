function [pos_W] = transformFromOtoW(windDirection_rad,pos_O)
%TRANSFORMFROMOTOW Summary of this function goes here
%   Detailed explanation goes here


M_WO = [cos(windDirection_rad), sin(windDirection_rad), 0;
        sin(windDirection_rad), -cos(windDirection_rad), 0;
        0, 0, -1];

pos_W = M_WO*pos_O;

end
