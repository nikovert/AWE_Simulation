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

function [ sol,p_C_W] = getTargetOnBoothNewton2(Lem, p_kite_W, c0, l_tether, direction)
%GETVTONLEMNISCATE Summary of this function goes here
%   Detailed explanation goes here
p_kite_W = p_kite_W/norm(p_kite_W);

[sol, ~] = doNewtonIterationBooth2(c0,Lem, p_kite_W, direction);

sol = sol + direction*0*pi/180;

long = Lem.b * sin(sol) ./( 1+(Lem.a/Lem.b*cos(sol)).^2 );
lat =   Lem.a * sin(sol).*cos(sol) ./ ( 1+(Lem.a/Lem.b*cos(sol)).^2 ) ;
L = [long;lat];
p_C_W = [cos(L(1,:)).*cos(L(2,:));
    sin(L(1,:)).*cos(L(2,:));
    sin(L(2,:))]*l_tether;


end