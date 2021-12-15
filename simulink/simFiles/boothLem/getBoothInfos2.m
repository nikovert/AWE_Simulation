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


function [t,DtDs,L, dLds,q] = getBoothInfos2(s_old,Lbooth, direction)
%GETBOOTHINFOS Calculate the tangent at the solution and its derivative

a = Lbooth.a; 
b = Lbooth.b;
L = [ b * sin(s_old) ./( 1+(a/b*cos(s_old)).^2 );
    a * sin(s_old).*cos(s_old) ./ ( 1+(a/b*cos(s_old)).^2 ) ] ;

dLds = [ ( b^3*cos(s_old).*(2*a^2-a^2*cos(s_old).^2+b^2)./(a^2*cos(s_old).^2+b^2).^2 );
    ( (cos(s_old).^2*(a^3*b^2+2*a*b^4) - a*b^4)./(a^2*cos(s_old).^2+b^2).^2 )];  

d2Lds2 =  [ ( -( (a^4*b^3*sin(s_old).^5-b^3*sin(s_old).*(5*a^4+4*(a*b)^2-b^4)+b^3*sin(s_old).^3 *(4*a^4+6*(a*b)^2 ))./(b^2-a^2*(sin(s_old).^2 - 1 ) ).^3 ) );
            ( (2*a*b^2*cos(s_old).^3.*sin(s_old)*(a^4+2*(a*b)^2)-a*b^2*cos(s_old).*sin(s_old)*(3*(a*b)^2+2*b^4)*2)./(a^2*cos(s_old).^2 + b^2 ).^3 )];


s_lambda = sin( L(1,:)  ); 
s_phi = sin( L(2,:) );
c_lambda = cos( L(1,:) ); 
c_phi = cos( L(2,:) ); 

q = [c_lambda*c_phi; s_lambda*c_phi; s_phi];

dqdlambda = [-s_lambda*c_phi; c_lambda*c_phi; 0];
dqdphi = [-c_lambda*s_phi; -s_lambda*s_phi; c_phi];

t = direction*( dqdlambda * dLds(1) + dqdphi *dLds(2) ); 

dtdlambda = [-c_lambda*c_phi*dLds(1)+s_lambda*s_phi*dLds(2); 
             -s_lambda*c_phi*dLds(1)-s_phi*c_lambda*dLds(2);
            0];
        
dtdphi = [s_lambda*s_phi*dLds(1)-c_phi*c_lambda*dLds(2); 
         -c_lambda*s_phi*dLds(1)-c_phi*s_lambda*dLds(2); 
         -s_phi*dLds(2)];
     
dtds = dqdlambda * d2Lds2(1)  + dqdphi * d2Lds2(2); 
     
% The negative sign cancels out 
DtDs =  ( dtdlambda * dLds(1) + dtdphi * dLds(2) + dtds ); 



end

