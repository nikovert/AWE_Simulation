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
% :Adapted from: Sebastian Rapp (s.rapp@tudelft.nl)

function [s_new,exceedMaxIter, nonconverged_points_list] = doNewtonIterationBooth2(s_old,Lem, pos_W, direction)
    res = 1;
    cnt = 1; 
    maxIter = 100; 
    s_new = s_old;
    nonconverged_points = inf;
    cnt_prev = 0;
    while max(res, [], 'all') > 0.01*pi/180 && cnt < maxIter
        if cnt > 2
            nonconverged_points = sum(res(:) > 0.01);
        end
        
        [t,DtDs] = getBoothInfos2(s_old,Lem, direction);    

        deltaDs =  pos_W{1} .* t{1} + pos_W{2} .* t{2} + pos_W{3} .* t{3}; 

        deltaD2s = pos_W{1} .* DtDs{1} + pos_W{2} .* DtDs{2} + pos_W{3} .* DtDs{3}; 

        s_new = s_old - direction * deltaDs ./ deltaD2s; % Note: We have to take into account the sign of the tangent
        % for the Minim. the direction of flight does not matter, we changed
        % however the sign according to it, hence we have to adapt it here again.

        s_new = mod(s_new, 2*pi);

        res = abs( s_new - s_old ); 

        s_old = s_new;
        
        % move s_old if we are stuck
        if nonconverged_points - sum(res(:) > 0.01) == 0 && cnt > cnt_prev+50
            index = find(res > 0.01);
            s_old(index) = s_old(min(index+1, length(s_old)));
            nonconverged_points = nan;
            cnt_prev = cnt;
        end
        
        cnt = cnt + 1 ;

    end
    
    nonconverged_points_list = find(res > 0.01);
    
    if cnt >= maxIter || nonconverged_points > 0
        exceedMaxIter = 1; 
    else
        exceedMaxIter = 0;
    end
end