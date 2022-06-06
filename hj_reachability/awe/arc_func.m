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