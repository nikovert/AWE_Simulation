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
% :Adapted from Mo Chen (mochen@sfu.ca)

function x1 = updateState(obj, u, T, x0, d)
% x1 = updateState(obj, u, T, x0, d)
% Updates state based on control
%
% Inputs:   obj - current quardotor object
%           u   - control (defaults to previous control)
%           T   - duration to hold control
%           x0  - initial state (defaults to current state if set to [])
%           d   - disturbance (defaults to [])
%
% Outputs:  x1  - final state
%

% If no state is specified, use current state
if nargin < 4 || isempty(x0)
  x0 = obj.x;
end

% If time horizon is 0, return initial state
if T == 0
  x1 = x0;
  return
end

% Default disturbance
if nargin < 5
  d = [];
end

% Do nothing if control is empty
if isempty(u)
  x1 = x0;
  return;
end

% convert u to vector if needed
if iscell(u)
  u = cell2mat(u);
end

if iscell(d)
  u = cell2mat(d);
end

% Do nothing if control is not a number
if isnan(u)
  warning('u = NaN')
  x1 = x0;
  return;
end

if isnan(d)
  warning('d = NaN')
  x1 = x0;
  return;
end

% Make sure control input is valid
if ~isnumeric(u)
  error('Control must be numeric!')
end

if ~isnumeric(d)
  error('Disturbance must be numeric!')
end

% Convert control to column vector if needed
if ~iscolumn(u)
  u = u';
end

if ~iscolumn(d)
  d = d';
end

% Check whether there's disturbance (this is needed since not all vehicle
% classes have dynamics that can handle disturbance)
options = odeset('RelTol',1e-5,'AbsTol',1e-7, 'MaxStep', 0.1);
if isempty(d)
  [~, x] = ode113(@(t,x) obj.dynamics(t, x, u), [0 T], x0, options);
else
  [~, x] = ode113(@(t,x) obj.dynamics(t, x, u, d), [0 T], x0, options);
end

% Update the state, state history, control, and control history
x1 = x(end, :)';
obj.x = x1;
obj.u = u;

obj.xhist = cat(2, obj.xhist, x1);
obj.uhist = cat(2, obj.uhist, u);
end