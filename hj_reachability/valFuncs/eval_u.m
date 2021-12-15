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
% :Author: Mo Chen (mochen@sfu.ca)

function v = eval_u(gs, datas, xs, interp_method)
% v = eval_u(g, datas, x)
%   Computes the interpolated value of the value functions datas at the
%   states xs
%
% Inputs:
%   Option 1: Single grid, single value function, multiple states
%     gs    - a single grid structure
%     datas - a single matrix (look-up stable) representing the value
%             function
%     xs    - set of states; each row is a state
%
%   Option 2: Single grid, multiple value functions, single state
%     gs    - a single grid structure
%     datas - a cell structure of matrices representing the value function
%     xs    - a single state
%
%   Option 3: Multiple grids, value functions, and states. The number of
%             grids, value functions, and states must be equal under this
%             option
%     gs    - a cell structure of grid structures
%     datas - a cell structure of matrices representing value functions
%     xs    - a cell structure of states
%
%   interp_method:  The default is linear interpolation.  
%                   Available methods are:
%       'nearest' - nearest neighbor interpolation
%       'linear'  - linear interpolation
%       'spline'  - spline interpolation
%       'cubic'   - cubic interpolation as long as the data is uniformly
%                   spaced, otherwise the same as 'spline'
%       'makima'  - modified Akima cubic interpolation
%
 
if nargin < 4
  interp_method = 'linear';
end
 
if isstruct(gs) && isnumeric(datas) && ismatrix(xs)
  % Option 1
  v = eval_u_single(gs, datas, xs, interp_method);
  
elseif isstruct(gs) && iscell(datas) && isvector(xs)
  % Option 2
  v = zeros(length(datas), 1);
  for i = 1:length(datas)
    v(i) = eval_u_single(gs, datas{i}, xs, interp_method);
  end
  
elseif iscell(gs) && iscell(datas) && iscell(xs)
  % Option 3
  v = zeros(length(gs), 1);
  for i = 1:length(gs)
    v(i) = eval_u_single(gs{i}, datas{i}, xs{i}, interp_method);
  end
  
else
  error('Unrecognized combination of input data types!')
end
end
 
function v = eval_u_single(g, data, x, interp_method)
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
%
% Mo Chen, 2015-10-15
% Updated 2016-05-18
 
% If the number of columns does not match the grid dimensions, try taking
% transpose
if size(x, 2) ~= g.dim
  x = x';
end
 
[g, data] = augmentPeriodicData(g, data);
 
%% Dealing with periodicity
for i = 1:g.dim
  if isfield(g, 'bdry') && isequal(g.bdry{i}, @addGhostPeriodic)
    % Map input points to within grid bounds
    period = max(g.vs{i}) - min(g.vs{i});
    
    i_above_bounds = x(:,i) > max(g.vs{i});
    while any(i_above_bounds)
      x(i_above_bounds, i) = x(i_above_bounds, i) - period;
      i_above_bounds = x(:,i) > max(g.vs{i});
    end
    
    i_below_bounds = x(:,i) < min(g.vs{i});
    while any(i_below_bounds)
      x(i_below_bounds, i) = x(i_below_bounds, i) + period;
      i_below_bounds = x(:,i) < min(g.vs{i});
    end
  end
end
 
%% Interpolate
% Input checking
x = checkInterpInput(g, x);
 
% eg. v = interpn(g.vs{1}, g.vs{2}, data, x(:,1), x(:,2), interp_method)
interpn_argin_x = cell(g.dim, 1);
for i = 1:g.dim
  interpn_argin_x{i} = x(:,i);
end
 
v = interpn(g.vs{:}, data, interpn_argin_x{:}, interp_method);
 
end
