function [p_unc] = addAdditiveAndMultiplUnc(p,scale, bias)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

p_unc = p*scale + bias; 

end

