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

function datatOut = regrid(dataIN, gridIn, gridOut, interp_method, tau)
% projects data onto a new grid
    if nargin < 4
        interp_method = 'linear';
    end
    
    assert(gridOut.dim == gridIn.dim)
    
    
    if nargin > 4
        gridIn.dim = gridIn.dim + 1;
        gridIn.vs{end+1} = tau;
        gridIn.shape = [gridIn.shape length(tau)];
        gridIn.bdry{end+1} = @addGhostExtrapolate;
        
        gridOut.dim = gridOut.dim + 1;
        gridOut.vs{end+1} = tau;
        gridOut.shape = [gridOut.shape length(tau)];
        points = zeros(prod(gridOut.shape), gridOut.dim);
        for i = 1:gridIn.dim
            x_out = full(gridOut.vs{i});
            s_out = ones(1,gridOut.dim); 
            s_out(i) = numel(x_out);
            x_out = reshape(x_out,s_out);
            s_out = gridOut.shape; 
            s_out(i) = 1;
            gridOut.xs{i} = repmat(x_out,s_out);
            clear x_out
            
            points(:, i) = gridOut.xs{i}(:);
            gridOut.xs{i} = [];
        end
        
    else
        points = zeros(prod(gridOut.shape), gridOut.dim);
        for i = 1:gridOut.dim
            points(:, i) = gridOut.xs{i}(:);
        end
    end
    v = eval_u(gridIn, dataIN, points, interp_method);
    clear dataIn points
    
    if isnan(v)
        warning('boundaries do not match, regridding unsuccessfull')
    end
    datatOut = reshape(v, gridOut.shape);
end