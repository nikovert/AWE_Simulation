function u = get_path_u(x, alpha_options, mu_options, alpha_max, grid_min, grid_max, dx, mu_max, I_table, alpha_min, mu_min)
    boundary_min = grid_min' > x;
    boundary_max = grid_max' < x;
    if any(boundary_min)
        warning('out of bounds')
        x(boundary_min) = grid_min(boundary_min);
    elseif any(boundary_max)
        warning('out of bounds')
        x(boundary_max) = grid_max(boundary_max);
    end

    I = eval_u(grid_min, grid_max, dx, I_table, x, 'nearest');
    [I1,I2] = ind2sub([alpha_options, mu_options], I);
    
    alpha_a = alpha_min + (I1-1)/(alpha_options-1) * (alpha_max - alpha_min);
    mu_a    = mu_min    + (I2-1)/(mu_options   -1) * (mu_max    - mu_min);
    u = [alpha_a, mu_a];
end

function v = eval_u(grid_min, grid_max, dx, datas, xs, interp_method)
% Only supporting Single grid, single value function, multiple states
    if nargin < 4
      interp_method = 'linear';
    end
    v = eval_u_single(grid_min, grid_max, dx, datas, xs, interp_method);
end

function v = eval_u_single(grid_min, grid_max, dx, data, x, interp_method)
% v = eval_u_single(g, data, x)
%   Computes the interpolated value of a value function data at state x
%
% Inputs:
%   g       - grid
%   data    - implicit function describing the set
%   x       - points to check; each row is a point
%
% OUTPUT
%   v:  value at points x


%% Dealing with periodicity
% Not implementing now, but should be handled

%% Interpolate
% Input checking
%x = checkInterpInput(g, x);

% Hardcoded for 7 states for now
VS_1 = (grid_min(1) : dx(1) : grid_max(1))';
VS_2 = (grid_min(2) : dx(2) : grid_max(2))';
VS_3 = (grid_min(3) : dx(3) : grid_max(3))';
VS_4 = (grid_min(4) : dx(4) : grid_max(4))';
VS_5 = (grid_min(5) : dx(5) : grid_max(5))';
VS_6 = (grid_min(6) : dx(6) : grid_max(6))';
VS_7 = (grid_min(7) : dx(7) : grid_max(7))';

v = interpn(VS_1, VS_2, VS_3, VS_4, VS_5, VS_6, VS_7, data, ...
            x(1), x(2), x(3), x(4), x(5), x(6), x(7), interp_method);

end