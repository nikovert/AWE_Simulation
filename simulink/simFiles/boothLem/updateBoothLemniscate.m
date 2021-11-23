function [Lem] = updateBoothLemniscate(l_tether,Lbooth)


%Lem.a = Lbooth.a/(l_tether*cos(Lbooth.phi0)); 
%Lem.b = Lbooth.b/(l_tether*cos(Lbooth.phi0)); 
Lem.a = Lbooth.a/(l_tether); 
Lem.b = Lbooth.b/(l_tether); 
Lem.phi0 = Lbooth.phi0;
