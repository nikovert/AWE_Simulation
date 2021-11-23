function [ sol,p_C_W] = getTargetOnBoothNewton2(Lem, p_kite_W, c0, l_tether, direction)
%GETVTONLEMNISCATE Summary of this function goes here
%   Detailed explanation goes here
p_kite_W = p_kite_W/norm(p_kite_W);

[sol, exceedMaxIter] = doNewtonIterationBooth2(c0,Lem, p_kite_W, direction);

sol = sol + direction*0*pi/180;

%L = [LemPs.Alambda * sin(LemPs.blambda*sol');
 %   LemPs.Aphi    * sin(LemPs.bphi*sol') + LemPs.phi0];
long = Lem.b * sin(sol) ./( 1+(Lem.a/Lem.b*cos(sol)).^2 );
lat =   Lem.a * sin(sol).*cos(sol) ./ ( 1+(Lem.a/Lem.b*cos(sol)).^2 ) ;
L = [long;lat];
p_C_W = [cos(L(1,:)).*cos(L(2,:));
    sin(L(1,:)).*cos(L(2,:));
    sin(L(2,:))]*l_tether;


end