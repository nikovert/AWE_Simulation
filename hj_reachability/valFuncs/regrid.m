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