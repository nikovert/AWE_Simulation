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

function [gOut, dataOut] = proj(g, data, dimsToRemove, xs, NOut, process)
% [gOut, dataOut] = proj(g, data, dims, xs, NOut)
%   Projects data corresponding to the grid g in g.dim dimensions, removing
%   dimensions specified in dims. If a point is specified, a slice of the
%   full-dimensional data at the point xs is taken.
%
% Inputs:
%   g            - grid corresponding to input data
%   data         - input data
%   dimsToRemove - vector of length g.dim specifying dimensions to project
%                 For example, if g.dim = 4, then dims = [0 0 1 1] would
%                 project the last two dimensions
%   xs           - Type of projection (defaults to 'min')
%       'min':    takes the union across the projected dimensions
%       'max':    takes the intersection across the projected dimensions
%       a vector: takes a slice of the data at the point xs
%   NOut    - number of grid points in output grid (defaults to the same
%             number of grid points of the original grid in the unprojected
%             dimensions)
%   process            - specifies whether to call processGrid to generate
%                        grid points
%
% Outputs:
%   gOut    - grid corresponding to projected data
%   dataOut - projected data
%
% See proj_test.m

%% Input checking
g.dim = length(g.shape);
% if g.shape(end) == 1
%     g.dim = g.dim -1;
% end
if length(dimsToRemove) ~= g.dim
  error('Dimensions are inconsistent!')
end
remove_dims = find(g.shape==1);
g.dim = g.dim - length(remove_dims);
dimsToRemove(remove_dims) = [];
g.xs(remove_dims) = [];
g.shape(remove_dims) = [];
data = reshape(data, g.shape);
if nnz(~dimsToRemove) == g.dim
  gOut = g;
  dataOut = data;
  warning('Input and output dimensions are the same!')
  return
end

% By default, do a projection
if nargin < 4
  xs = 'min';
end

% If a slice is requested, make sure the specified point has the correct
% dimension
if isnumeric(xs) && length(xs) ~= nnz(dimsToRemove)
  error('Dimension of xs and dims do not match!')
end

if nargin < 5
  NOut = g.N(~dimsToRemove);
end

if nargin < 6
  process = true;
end

dataDims = numDims(data);
if ~isempty(data) && ~(dataDims == g.dim || dataDims == g.dim+1) && ~iscell(data)
  error('Inconsistent input data dimensions!')
end

%% Project data
if dataDims == g.dim
  [gOut, dataOut] = projSingle(g, data, dimsToRemove, xs, NOut, process);
  
else % dataDims == g.dim + 1
  % Project grid
  gOut = projSingle(g, [], dimsToRemove, xs, NOut, process);
  
  % Project data
  if iscell(data)
    numTimeSteps = length(data);
  else
    numTimeSteps = size(data, dataDims);
    colonsIn = repmat({':'}, 1, g.dim);
  end
  
  dataOut = zeros([NOut' numTimeSteps]);
  colonsOut = repmat({':'}, 1, gOut.dim);
  
  for i = 1:numTimeSteps
    if iscell(data)
    [~, dataOut(colonsOut{:},i)] = ...
      projSingle(g, data{i}, dimsToRemove, xs, NOut, process);   
    else
    [~, dataOut(colonsOut{:},i)] = ...
      projSingle(g, data(colonsIn{:},i), dimsToRemove, xs, NOut, process);      
    end
  end
end

end
function [gOut, dataOut] = projSingle(g, data, dims, xs, NOut, process)
% [gOut, dataOut] = proj(g, data, dims, xs, NOut)
%   Projects data corresponding to the grid g in g.dim dimensions, removing
%   dimensions specified in dims. If a point is specified, a slice of the
%   full-dimensional data at the point xs is taken.
%
% Inputs:
%   g       - grid corresponding to input data
%   data    - input data
%   dims    - vector of length g.dim specifying dimensions to project
%                 For example, if g.dim = 4, then dims = [0 0 1 1] would
%                 project the last two dimensions
%   xs      - Type of projection (defaults to 'min')
%       'min':    takes the union across the projected dimensions
%       'max':    takes the intersection across the projected dimensions
%       a vector: takes a slice of the data at the point xs
%   NOut    - number of grid points in output grid (defaults to the same
%             number of grid points of the original grid in the unprojected
%             dimensions)
%   process            - specifies whether to call processGrid to generate
%                        grid points
%
% Outputs:
%   gOut    - grid corresponding to projected data
%   dataOut - projected data
%
% See proj_test.m

%% Create ouptut grid by keeping dimensions that we are not collapsing
if isempty(g)
  if ~ischar(xs) || (~strcmp(xs, 'max') && ~strcmp(xs, 'min'))
    error('Must perform min or max projection when not specifying grid!')
  end
  
else
  dims = logical(dims);
  gOut.dim = nnz(~dims);
  gOut.min = g.min(~dims);
  gOut.max = g.max(~dims);
  gOut.bdry = g.bdry(~dims);
  
  if numel(NOut) == 1
    gOut.N = NOut*ones(gOut.dim, 1);
  else
    gOut.N = NOut;
  end
  
  % Process the grid to populate the remaining fields if necessary
  if process
    gOut = processGrid(gOut);
  end
  
  % Only compute the grid if value function is not requested
  if nargout < 2
    return
  end
end

%% 'min' or 'max'
if ischar(xs)
  dimsToProj = find(dims);
  
  for i = length(dimsToProj):-1:1
    if strcmp(xs,'min')
      data = squeeze(min(data, [], dimsToProj(i)));
    elseif strcmp(xs,'max')
      data = squeeze(max(data, [], dimsToProj(i)));
    else
      error('xs must be a vector, ''min'', or ''max''!')
    end
  end
  
  dataOut = data;
  return
end

%% Take a slice
% Preprocess periodic dimensions
[g, data] = augmentPeriodicData(g, data);

% temp = interpn(g.vs{1}, g.vs{2}, g.vs{3}, g.vs{4}, data, g.vs{1}, xs(1), ...
%   g.vs{3}, xs(2));
eval_pt = cell(g.dim, 1);
xsi = 1;
for i = 1:g.dim
  if dims(i)
    % If this dimension is periodic, wrap the input point to the correct period
    if isfield(g, 'bdry') && isequal(g.bdry{i}, @addGhostPeriodic)
      period = max(g.vs{i}) - min(g.vs{i});
      
      while xs(xsi) > max(g.vs{i})
        xs(xsi) = xs(xsi) - period;
      end
      
      while xs(xsi) < min(g.vs{i})
        xs(xsi) = xs(xsi) + period;
      end
    end
    
    eval_pt{i} = xs(xsi);
    xsi = xsi + 1;
    
  else
    eval_pt{i} = g.vs{i};
  end
  
end
temp = interpn(g.vs{:}, data, eval_pt{:});

dataOut = squeeze(temp);

% interpn(g.vs{1}, g.vs{3}, dataOut, gOut.xs{1}, gOut.xs{2})
dataOut = interpn(g.vs{~dims}, dataOut, gOut.xs{:});
end