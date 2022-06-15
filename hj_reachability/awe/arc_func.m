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
% :Revision: 15-June-2021
% :Author: Nikolaus Vertovec (nikolaus.vertovec@eng.ox.ac.uk)

function [points] = arc_func(PoI, plane, extraArgs)
    if nargin < 3
        extraArgs = [];
    end
    
    if isfield(extraArgs, 'scale')
        scale = extraArgs.scale;
    else
        scale = 10;
    end

    if isfield(extraArgs, 'nrpoints')
        nrpoints = extraArgs.nrpoints;
    else
        nrpoints = 500;
    end

    if isfield(extraArgs, 'M')
        M = extraArgs.M;
        if iscell(M)
            M = cell2mat(M);
        end
    else
        M = eye(3);
    end
    if iscell(PoI)
        PoI = cell2mat(PoI);
    end
    switch plane
        case {'xy','yx'}
            M = @(t) [cos(t), -sin(t), 0; ...
                     sin(t), cos(t), 0; ...
                     0, 0, 1];
        case {'xz','zx'}
            M = @(t) [cos(t), 0, sin(t); ...
                     0, 1, 0; ...
                     -sin(t), 0, cos(t)];
        case {'yz','zy'}
            M = @(t) [1, 0, 0; ...
                     0, cos(t), -sin(t); ...
                     0, sin(t), cos(t)];
        otherwise
            error('unknown plane')
    end
    points = zeros(3, nrpoints);
    for i = 1:nrpoints
        points(:,i) = PoI + M(i * 2 *pi/nrpoints) * M * [scale; 0; 0];
    end
    
end
