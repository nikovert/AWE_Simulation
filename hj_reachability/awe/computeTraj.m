function [traj, traj_tau] = computeTraj(g, data, tau, dynSys, extraArgs)
% [traj, traj_tau] = computeOptTraj(g, data, tau, dynSys, extraArgs)
%   Computes the optimal trajectories given the optimal value function
%   represented by (g, data), associated time stamps tau, dynamics given in
%   dynSys.
%
% Inputs:
%   g, data - grid and value function
%   tau     - time stamp (must be the same length as size of last dimension of
%                         data)
%   dynSys  - dynamical system object for which the optimal path is to be
%             computed
%   extraArgs
%     .uMode        - specifies whether the control u aims to minimize or
%                     maximize the value function
%     .visualize    - set to true to visualize results
%     .fig_num:   List if you want to plot on a specific figure number
%     .projDim      - set the dimensions that should be projected away when
%                     visualizing
%     .fig_filename - specifies the file name for saving the visualizations
%     .pdDims

if nargin < 5
  extraArgs = [];
end

% Default parameters
vis = false;
subSamples = 10;

if isfield(extraArgs, 'uMode')
  uMode = extraArgs.uMode;
else
  uMode = 'min';
end

if isfield(extraArgs, 'dMode')
  dMode = extraArgs.dMode;
else
  dMode = 'max';
end

if isfield(extraArgs, 'derivFunc')
  derivFunc = extraArgs.derivFunc;
else
  derivFunc = @upwindFirstWENO5;
end

% Visualization
if isfield(extraArgs, 'visualize') && extraArgs.visualize
  vis = extraArgs.visualize;
  
  showDims = find(extraArgs.projDim);
  hideDims = ~extraArgs.projDim;
  
  if isfield(extraArgs,'fig_num')
    fig_num = extraArgs.fig_num;
  else
    fig_num = 2;
  end
  set_fig = figure(fig_num);
  traj_fig = figure(fig_num+1);
  force_fig = figure(fig_num+2);
end

if isfield(extraArgs, 'subSamples')
  subSamples = extraArgs.subSamples;
end

if isfield(extraArgs, 'interpolationMethod')
  interpolationMethod = extraArgs.interpolationMethod;
else
  interpolationMethod = 'linear';
end

if isfield(extraArgs, 'continuous')
  continuous = extraArgs.continuous;
else
  continuous = false;
end

clns = repmat({':'}, 1, g.dim);

if any(diff(tau) < 0)
  error('Time stamps must be in ascending order!')
end

% Time parameters
tauLength = length(tau);
dtSmall = (tau(2) - tau(1))/subSamples;
% maxIter = 1.25*tauLength;

% Initialize trajectory
traj = nan(g.dim, tauLength);
traj(:,1) = dynSys.x;
tEarliest = 1;

upper = tauLength;
lower = tEarliest;
iter = 1;

g.pdDims = [];
for i=1:g.dim
    if isequal(g.bdry{i}, @addGhostPeriodic)
        g.pdDims = [g.pdDims i];
    end
end

while max(tEarliest,iter) <= tauLength || continuous
  if continuous
      tauLength = tauLength + 1;
  end
  % Determine the earliest time that the current state is in the reachable set
  % Binary search
  
  tEarliest = find_earliest_BRS_ind(g, data, dynSys.x, upper, lower);
  %tEarliest = iter;
  % BRS at current time
  BRS_at_t = data(clns{:}, tEarliest);
  
  % Visualize BRS corresponding to current trajectory point
  if vis
    set(0, 'CurrentFigure', set_fig)
    plot(traj(showDims(1), iter), traj(showDims(2), iter), 'kh')
    hold on
    [g2D, data2D] = proj(g, BRS_at_t, hideDims, traj(hideDims,iter));
    if any(isnan(data2D))
       warning('Projection unsuccesful') 
    end
    visSetIm(g2D, data2D);
    tStr = sprintf('t = %.3f; tEarliest = %.3f', (iter-1)*dtSmall, tau(tEarliest));
    title(tStr)
    
    if isfield(extraArgs, 'fig_filename')
      export_fig(sprintf('%s%d', extraArgs.fig_filename, iter), '-png')
    end

    %hold off
  end
  
  if tEarliest == tauLength
    % Trajectory has entered the target
    break
  end
  
  % Update trajectory
  Deriv = computeGradients(g, BRS_at_t, true(g.dim, 1), derivFunc);
  for j = 1:subSamples
    deriv = eval_u(g, Deriv, dynSys.x, interpolationMethod);
    if any(isnan(deriv))
        deriv = eval_u(g,Deriv, dynSys.x, 'spline');
    end
    u = dynSys.optCtrl(tau(tEarliest), dynSys.x, deriv, uMode);
    d = dynSys.optDstb(tau(tEarliest), dynSys.x, deriv, dMode);
    %d = (2*rand-1)*dynSys.max_tether_diff_dot;
    dynSys.updateState(u, dtSmall, dynSys.x, d);
    
    boundary_min = g.min >= dynSys.x;
    boundary_max = g.max <= dynSys.x;
    boundary_min(g.pdDims) = [];
    boundary_max(g.pdDims) = [];
    if any(boundary_min) || any(boundary_max)
        warning('out of bounds')
        iter = iter + 1;
        traj(:,iter) = dynSys.x;
        traj(:,iter+1:end) = [];
        traj_tau = tau(1:min(iter-1, length(tau)));
        return
    end
  end
  if vis
    extraArgs.u = u;
    extraArgs.fig_handle = traj_fig;
    extraArgs.visualizePath = true;
    dynSys.getTargetdistance(dynSys.x, extraArgs);
    
    [~, ~, max_F_tether] = dynSys.dynamics(0, dynSys.x, u);
    Ft_ground_norm = norm_cellVec(max_F_tether)/1e+03;
    set(0, 'CurrentFigure', force_fig)
    hold on
    plot(iter, Ft_ground_norm, 'kh')
    if Ft_ground_norm > dynSys.F_T_max
        warning('Tether rupture')
    end
  end
  % Record new point on nominal trajectory
  iter = iter + 1;
  traj(:,iter) = dynSys.x;
      
end

% Delete unused indices
traj(:,iter+1:end) = [];
traj_tau = tau(1:min(length(tau),iter-1));
end

function tEarliest = find_earliest_BRS_ind(g, data, x, upper, lower)
% tEarliest = find_earliest_BRS_ind(g, data, x, upper, lower)
%     Determine the earliest time that the current state is in the reachable set
% 
% Inputs:
%     g, data - grid and value function representing reachable set
%     x       - state of interest
%     upper, lower - upper and lower indices of the search range
%
% Output:
%     tEarliest - earliest time index that x is in the reachable set

if nargin < 4
  upper = size(data, g.dim+1);
  lower = 1;
end

% Binary search
clns = repmat({':'}, 1, g.dim);
small = 1e-6;

while upper > lower
  tEarliest = ceil((upper + lower)/2);
  valueAtX = eval_u(g, data(clns{:}, tEarliest), x);
  
  if valueAtX < small
    % point is in reachable set; eliminate all lower indices
    lower = tEarliest;
  else
    % too late
    upper = tEarliest - 1;
  end
end

tEarliest = upper;
end