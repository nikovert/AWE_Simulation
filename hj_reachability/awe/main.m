 %% Readme
% This script is intended to calculate the optimal path for tracking a
% figure 8 using a 3 DOF AWE model
%
% Last updated 23 Nov 2021
%%
% clear;
% close all;
% clc;

%% Setup

%% Add relevant files to path
addpath('../')
addToolbox;
    
filename = 'main_run.mat';

%% Settings
visualize_contour = true;
%% Setup sys

% Normalisation parameter
h0 = 1;
v0 = 1;
a0 = 1;

s       = 4*pi/3/a0;
sigma   = 0;
h_tau   = 250/h0;
Va      = 31/v0;
chi_a   = -0.4470/a0;
gamma_a = 0.5205/a0;
tether_diff = 0.0003;
initialState = [s, sigma, h_tau, Va, chi_a, gamma_a, tether_diff]';

%% Grid 
%                long,    lat,   h_tau,     Va,     chi_a,    gamma_a,  tether_diff (not normalised)
N        = [       31;   7;       7;      9;        10;         11;     9];
%N        = [        3;   3;       3;      3;        3;         3;     3];
grid_min = [     0/a0; -45;  200/h0;  20/v0;    -pi/a0;   -pi/3/a0; -1e-3]; 
grid_max = [  2*pi/a0;  45;  600/h0;  40/v0;     pi/a0;    pi/3/a0;  7e-3];
pdDims   = [1 5];
process  = true;
% N = N(1:end-1);
% grid_min = grid_min(1:end-1);
% grid_max = grid_max(1:end-1);
grid = createGrid(grid_min, grid_max, N, pdDims, process);
    
disp(['using following grid: ', num2str(grid.shape)])
%% Setup dynSys
sys = AWE_3DOF(initialState(1:grid.dim));
sys.h0 = h0;
sys.v0 = v0;
sys.a0 = a0;
sys.v_w_O = {-15.8031/v0   -0.0001/v0   -0.0000/v0}';
sys.Ft_set = 1.652199337711693e+03;
sys.v_ro_set = 0.01; % don't normalise
sys.F_T_max = 1.850;
sys.skipTether = false; % skips the calculation of the tether force
sys.ignoreTetherDiff = false; % assumes that the tether diff never changes
sys.curve_direction = -1;
Lem.a = 120;
Lem.b = 200;
Lem.phi0 = 1;
sys.Lem = Lem;
%% Set up mask for set K
alpha   = 0.0;
mu_a    = 0.0;
u = {alpha, mu_a};
[~, ~, max_F_tether] = sys.dynamics(0, grid.xs, u);
Ft_ground_norm = norm_cellVec(max_F_tether)/1e+03;

mask_force = (Ft_ground_norm - sys.F_T_max)/max(abs(Ft_ground_norm(:) - sys.F_T_max));

mask_bndry = -inf(size(mask_force));
for i = 1 : grid.dim
  if any(i==pdDims)
      continue
  end
  mask_bndry = max(mask_bndry, grid.xs{i} - grid.max(i));
  mask_bndry = max(mask_bndry, grid.min(i) - grid.xs{i});
end

mask = shapeIntersection(mask_force,mask_bndry);
mask = mask_force;
assert(eval_u(grid, mask, initialState(1:grid.dim), 'linear') <= 0);

%% Define Target States
targetDistanceArgs.Lem = Lem;
targetDistanceArgs.distanceOnly = false;
targetDistanceArgs.headingOnly = true;
targetDistanceArgs.normalize = false;
targetDistanceArgs.direction = sys.curve_direction;
%targetDistanceArgs.h_tau = 250;
[distance, sol,p_C_W, p_kite_W] = sys.getTargetdistance(grid.xs, targetDistanceArgs);
distance_sort = sort(distance(:));
data0 = distance-distance_sort(floor(length(distance_sort)/10)); % Make sure 10 percent of points are in the starting set

%% Clean up unused large sets
clear distance_sort distance mask_force mask_bndry Ft_ground_norm max_F_tether
%% Sanity check
if false
    long  = linspace(grid.min(1), grid.max(1), 60);
    lat   = linspace(grid.min(2), grid.max(2), 30);
    h_tau = grid.max(3) * ones(size(long));
    f = nan(length(long), length(lat));

    cmap = jet(100);
    colormap(cmap); colorbar;

    for long_index = 1:length(long)
        for lat_index = 1:length(lat)
            state = [long(long_index), lat(lat_index), h_tau(1), initialState(4:end)']';
            [d, sol,p_C_W, p_kite_W] = sys.getTargetdistance(state);
            f(long_index, lat_index) = d;
            d = max(0.05, min(d,1));
            c = cmap(ceil(100*d/1), :);
            scatter3(p_kite_W{1},p_kite_W{2},p_kite_W{3},...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor',c)
        end
        drawnow
    end
end

%% Pack problem parameters
% Put grid and dynamic systems into schemeData
schemeData.grid = grid;
schemeData.dynSys = sys;
schemeData.accuracy = 'medium'; %set accuracy
schemeData.uMode = 'min';
schemeData.dMode = 'max';
schemeData.dissType = 'local';
schemeData.hamFunc = @sys.Hamiltonian;
schemeData.partialFunc = @sys.partialFunc;

schemeData.tMode = 'backward';
%% Solver setup
if visualize_contour
    HJIextraArgs.visualize.valueFunction = 1;
    HJIextraArgs.visualize.valueSet = 1;

    HJIextraArgs.visualize.initialValueFunction = 1;
    HJIextraArgs.visualize.initialValueSet = 1;

    HJIextraArgs.visualize.figNum = 1; %set figure number
    HJIextraArgs.visualize.deleteLastPlot = true; %delete previous plot as you update

    % 2D slice
    HJIextraArgs.visualize.plotData.plotDims = [1, 1, zeros(1, grid.dim-2)]; %plot r, vt
    HJIextraArgs.visualize.plotData.projpt = initialState(3:grid.dim)';
else
    HJIextraArgs.visualize = false;
end

HJIextraArgs.makeVideo = false;
HJIextraArgs.videoFilename = 'value_function.mp4';
HJIextraArgs.frameRate = 10;

%% Extra Arguments
HJIextraArgs.saveFilename = strrep(filename, '.mat', '_tmp.mat');
HJIextraArgs.saveFrequency = 10;
HJIextraArgs.keepLast = true;
HJIextraArgs.obstacleFunction = -mask;

% Let the integrator know what function to call.
HJIextraArgs.lowMemory = false; % lowMemory is already implemented
HJIextraArgs.factorCFL = 0.8;
%% Set up time vector
t0   = 0;
tMax = 0.1;
dt   = 0.05; 
tau  = t0:dt:tMax;

%% Solve
[data, tau2, ~] = ...
    HJIPDE_solve(data0, tau, schemeData, 'set', HJIextraArgs);

% Save data
disp('Saving workspace')
save(filename, '-v7.3')

extraArgs.saveFile = true;
recast_to_new_grid = true;

if recast_to_new_grid
    N        = [       35;   9;       9;      9;        10;         11;     11]; 
    clear data0 mask sol
    extraArgs.newgrid = createGrid(grid_min, grid_max, N, pdDims);
end

if HJIextraArgs.keepLast
    [alpha_table, mu_table, I_table] = generate_lookup_table(grid, data, 1, sys, extraArgs);
    return
else
    [alpha_table, mu_table, I_table] = generate_lookup_table(grid, flip(data,length(grid.N)+1), 1, sys, extraArgs);
end
%% Trajectory
TrajextraArgs.visualize = true; %show plot
TrajextraArgs.fig_num = 2; %figure number
if grid.dim == 6
    dataf = data(:,:,:,:,:,:,end);
else
    dataf = data(:,:,:,:,:,:,:,end);
end
ind     = find(dataf< 0);       % All points in the BRT
index   = find(data0(ind)>0);   % index of int that is not in initial state
if isempty(index)
    return
end
I = cell(grid.dim,1);
if grid.dim == 6
    [I{1},I{2},I{3},I{4},I{5},I{6}]      = ind2sub(N,ind(index)); % convert index
else
    [I{1},I{2},I{3},I{4},I{5},I{6},I{7}] = ind2sub(N,ind(index)); % convert index
end
rm_ind = [];
for i=1:grid.dim % remove points on the boundary
    rm_ind = [rm_ind; find(I{i} == 1); find(I{i} == N(i))];
end
index(rm_ind) = [];
if grid.dim == 6
    [I{1},I{2},I{3},I{4},I{5},I{6}]      = ind2sub(N,ind(index)); % convert index
else
    [I{1},I{2},I{3},I{4},I{5},I{6},I{7}] = ind2sub(N,ind(index)); % convert index
end
for i=1:grid.dim
    sys.x(i) = grid.vs{i}(I{i}(ceil(2*length(index)/4)));
end

assert(eval_u(grid, dataf, sys.x, 'linear') < 0)
sys.xhist = [];
     
%we want to see the first two dimensions (x and y)
TrajextraArgs.projDim = HJIextraArgs.visualize.plotData.plotDims; 
TrajextraArgs.continuous = true;
TrajextraArgs.Lem = targetDistanceArgs.Lem;
TrajextraArgs.distanceOnly = targetDistanceArgs.distanceOnly;
TrajextraArgs.normalize = false;
if isfield(targetDistanceArgs, 'h_tau')
    TrajextraArgs.h_tau = targetDistanceArgs.h_tau;
end
dataTraj = flip(data,length(grid.N)+1);
[traj, traj_tau] = computeTraj(grid, dataTraj, tau, sys, TrajextraArgs);


