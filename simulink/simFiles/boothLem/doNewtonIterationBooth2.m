function [s_new,exceedMaxIter] = doNewtonIterationBooth2(s_old,Lem, pos_W, direction)

pos_W = pos_W/norm(pos_W); 
res = 1;
cnt = 1; 
maxIter = 10; 
s_new = s_old;
while res > 0.1*pi/180 && cnt < maxIter

[t,DtDs,L, dLds,q] = getBoothInfos2(s_old,Lem, direction);    
    
deltaDs =  pos_W'*t ; 

deltaD2s = pos_W'*DtDs; 

s_new = s_old - direction * deltaDs/deltaD2s; % Note: We have to take into account the sign of the tangent
% for the Minim. the direction of flight does not matter, we changed
% however the sign according to it, hence we have to adapt it here again.

s_new = mod(s_new, 2*pi);

res = abs( s_new - s_old ); 

s_old = s_new;


cnt = cnt + 1 ;

end

if cnt > maxIter 
    exceedMaxIter = 1; 
else
     exceedMaxIter = 0;
end

