%% Path setup
addpath('../')
addToolbox;
%% Normalisation parameter
h0 = 1;
v0 = 1;
a0 = 1;
direction = -1;

s       = 0*pi/a0;
sigma   = 0;
h_tau   = 250/h0;
Va      = 30/v0;
chi_a   = 0.7381/a0;
gamma_a = 0.2148/a0;
tether_diff = 0;
initialState = [s, sigma, h_tau, Va, chi_a, gamma_a, tether_diff]';
%% Setup dynSys
sys = AWE_3DOF(initialState);
sys.h0 = h0;
sys.v0 = v0;
sys.a0 = a0;
sys.v_w_O = {-15.8031/v0   -0.0001/v0   -0.0000/v0}';
sys.v_w_O = {0 0 0}'; % Update this later to check
sys.Ft_set = 1.652199337711693e+03;
sys.v_ro_set = 0.01; % don't normalise
sys.F_T_max = 1.6649;
sys.curve_direction = direction;

%% Define Target States
% Check that points on the curve actually lead to a distance of zero (distanceOnly = 1)
Lem.a = 120;
Lem.b = 200;
Lem.phi0 = 1;
sol_spread = 0:0.1:2*pi;

targetDistanceArgs.Lem = Lem;
targetDistanceArgs.distanceOnly = true;
targetDistanceArgs.visualize = true;
targetDistanceArgs.s_old = sol_spread;
targetDistanceArgs.normalize = false;
targetDistanceArgs.direction    = direction;

M_WP = {cos(Lem.phi0),0, -sin(Lem.phi0);0, 1, 0; sin(Lem.phi0),0, cos(Lem.phi0)};
Lem.a = Lem.a / (h0*h_tau);
Lem.b = Lem.b / (h0*h_tau);

long_path = Lem.b .* sin(sol_spread) ./ ( 1+(Lem.a./Lem.b .* cos(sol_spread)).^2 );
lat_path  = Lem.a .* sin(sol_spread) .* cos(sol_spread) ./ ( 1+(Lem.a./Lem.b .* cos(sol_spread)).^2 );
[X,Y,Z] = sph2cart(long_path,lat_path,h_tau);
points = cell2mat(M_WP) * [X;Y;Z];
[long, lat, r] = cart2sph(points(1,:), points(2,:), points(3,:)); 
    
State = cell(7,1);
State{1} = long/sys.a0;
State{2} = lat/sys.a0;
State{3} = h_tau;
State{4} = Va;
State{5} = chi_a;
State{6} = gamma_a;
State{7} = tether_diff;
targetDistanceArgs.LongLatState = true;
[distance, sol,p_C_W, p_kite_W] = sys.getTargetdistance(State, targetDistanceArgs);

if targetDistanceArgs.distanceOnly && any(distance > 1e-4)
    warning('error finding best point on curve')
end
[long_reverse, lat_reverse] = getLongLat(sol, distance, sys.h0 * h_tau);
if any(abs(long_reverse - long) > 1e-4) || any(abs(lat_reverse - lat) > 1e-4)
   warning('long_reverse and long as well as lat_reverse and lat should match') 
end
%% Single test
zeroState = [0.9058/a0; 0.5761/a0; 250/h0; 30/v0; -2.4578/a0; 0.0289/a0; 0];
targetDistanceArgs.LongLatState = true;
targetDistanceArgs.distanceOnly = false;
targetDistanceArgs.headingOnly = false;
[distance, sol,p_C_W, p_kite_W] = sys.getTargetdistance(zeroState, targetDistanceArgs);
if distance > 0.01
   warning('zeroState should lead to distance=0') 
end
%% Plot Spread
% Create a mesh of all the different long/lat values (including points on the curve)
% The points all the points on the grid should at some point be a target
[LONG, LAT] = meshgrid(linspace(min(long)-0.2, max(long)+0.2, 50),linspace(min(lat)-0.2, max(lat)+0.2, 50));
State{1} = LONG(:)'/sys.a0;
State{2} = LAT(:)'/sys.a0;

targetDistanceArgs.distanceOnly = true;
targetDistanceArgs.visualize = false;
targetDistanceArgs.LongLatState = true;
if isfield(targetDistanceArgs,'s_old')
    targetDistanceArgs = rmfield(targetDistanceArgs,'s_old');
end

[distance_total, sol_total,p_C_W_total, p_kite_W_total] = sys.getTargetdistance(State, targetDistanceArgs);
%surf(LONG, LAT, reshape(distance_total, size(LONG)))

figure2 = figure(2);

% Create axes
axes1 = axes('Parent',figure2);
hold(axes1,'on');

% Create contour
contourf(LONG, LAT, reshape(distance_total, size(LONG)));

% Create ylabel
ylabel('Latitude');

% Create xlabel
xlabel('Longitude');

box(axes1,'on');
axis(axes1,'tight');
hold(axes1,'off');
% Set the remaining axes properties
set(axes1,'BoxStyle','full','FontSize',20,'Layer','top');
% Create colorbar
colorbar(axes1);


figure(3)
% [X,Y,Z] = sph2cart(LONG(:),LAT(:),h_tau);
% mesh(reshape(X, size(LONG)), reshape(Y, size(LONG)), reshape(distance, size(LONG)))
surf(LONG, LAT, reshape(sol_total, size(LONG)))

%% Test individual points
% make sure that results of the previous subsection are the same when
% testing all the points individually
figure(4)
for i = 1:length(LONG(:))
    plot3(p_C_W_total{1}(i), p_C_W_total{2}(i), p_C_W_total{3}(i), 'b+'); hold on
    plot3(p_kite_W_total{1}(i), p_kite_W_total{2}(i), p_kite_W_total{3}(i), 'go');
    
    State{1} = LONG(i)'/sys.a0;
    State{2} = LAT(i)'/sys.a0;
    [distance_i, sol_i,p_C_W_i, p_kite_W_i] = sys.getTargetdistance(State, targetDistanceArgs);
    if targetDistanceArgs.distanceOnly && abs(distance_i - distance_total(i)) > 0.001
       warning('distance missmatch') 
    end
end

%% Test grid used in main
% test that the distance plot is equally smooth when using the grid as
% defined in main.m

%                 sol,  sigma,   h_tau,     Va,    chi_a,    gamma_a,  l_s_diff
N        = [       61;  3;       9;      9;        10;         11;     7]; 
grid_min = [     0/a0; -5;  240/h0;  25/v0;    -pi/a0;   -pi/3/a0; -1e-3]; 
grid_max = [  2*pi/a0;  5;  260/h0;  35/v0;     pi/a0;    pi/3/a0;  7e-3];
test_full_grid = true;
if test_full_grid
    pdDims   = [1 5];
    process  = true;
    grid = createGrid(grid_min(1:2), grid_max(1:2), N(1:2), pdDims(1), process);
    State{1} = grid.xs{1};
    State{2} = grid.xs{2};
    State{3} = h_tau;
    State{4} = Va;
    State{5} = chi_a;
    State{6} = gamma_a;
    State{7} = tether_diff;
    getLongLatArgs.visualize_steps = false;
    [long_grid, lat_grid] = getLongLat(sys.a0 .* grid.xs{1}, grid.xs{2}, sys.h0 * h_tau, getLongLatArgs);
else
    sol     = linspace(grid_min(1), grid_max(1), N(1))';
    sigma   = linspace(grid_min(2), grid_max(2), N(2))';
    [SOL, SIGMA] = meshgrid(sol,sigma);
    getLongLatArgs.visualize_steps = true;
    [long_grid, lat_grid] = getLongLat(sys.a0 .* SOL, SIGMA, sys.h0 * h_tau, getLongLatArgs);

    State{1} = SOL;
    State{2} = SIGMA;
    State{3} = h_tau;
    State{4} = Va;
    State{5} = chi_a;
    State{6} = gamma_a;
    State{7} = tether_diff;
end

targetDistanceArgs.visualize = false;
targetDistanceArgs.distanceOnly = false;
targetDistanceArgs.headingOnly = true;
targetDistanceArgs.LongLatState = false;
targetDistanceArgs.normalize = false;
if isfield(targetDistanceArgs,'s_old')
    targetDistanceArgs = rmfield(targetDistanceArgs,'s_old');
end
test_full_grid = false;
if test_full_grid
    grid = createGrid(grid_min(1:6), grid_max(1:6), N(1:6), pdDims, process);
    [distance_grid, sol_grid,p_C_W_grid, p_kite_W_grid] = sys.getTargetdistance(grid.xs, targetDistanceArgs);
    [gOut, dataOut] = proj(grid, distance_grid, ~[1, 1, zeros(1, grid.dim-2)], initialState(3:grid.dim)');
    figure(6)
    plot3(long_grid, lat_grid, zeros(size(State{1})), 'kh'); hold on
    h = surf(long_grid, lat_grid, reshape(dataOut, size(State{1})));
    h.EdgeColor = 'none'; 
    h.FaceLighting = 'phong';

    figure(7)
    plot3(long_grid, lat_grid, zeros(size(State{1})), 'kh'); hold on
    h = surf(gOut.xs{1}, gOut.xs{2}, reshape(dataOut, size(gOut.xs{1})));
    h.EdgeColor = 'none'; 
    h.FaceLighting = 'phong';
else
    figure(5)
    targetDistanceArgs.visualize = true;
    [distance_total, sol_total,p_C_W_total, p_kite_W_total] = sys.getTargetdistance(State, targetDistanceArgs);
    figure(6)
    plot3(long_grid, lat_grid, zeros(size(State{1})), 'kh'); hold on
    h = surf(long_grid, lat_grid, reshape(distance_total, size(State{1})));
    h.EdgeColor = 'none'; 
    h.FaceLighting = 'phong';

    figure(7)
    plot3(long_grid, lat_grid, zeros(size(State{1})), 'kh'); hold on
    h = surf(State{1}/sys.a0, State{2}, reshape(distance_total, size(State{1})));
    h.EdgeColor = 'none'; 
    h.FaceLighting = 'phong';
end

%% Test heading distance
%                 sol,  sigma,   h_tau,     Va,    chi_a,    gamma_a,  l_s_diff
N        = [       61;  3;       9;      9;        10;         11;     7]; 
grid_min = [     0/a0; -5;  240/h0;  25/v0;    -pi/a0;   -pi/3/a0; -1e-3]; 
grid_max = [  2*pi/a0;  5;  260/h0;  35/v0;     pi/a0;    pi/3/a0;  7e-3];
pdDims   = [1 5];
process  = true;

targetDistanceArgs.visualize = false;
targetDistanceArgs.distanceOnly = false;
targetDistanceArgs.headingOnly = false;
targetDistanceArgs.LongLatState = false;
targetDistanceArgs.normalize = false;
if isfield(targetDistanceArgs,'s_old')
    targetDistanceArgs = rmfield(targetDistanceArgs,'s_old');
end

grid = createGrid(grid_min(5:6), grid_max(5:6), N(5:6), 2, process);

State{1} = 0;
State{2} = 0;
State{3} = h_tau;
State{4} = Va;
State{5} = grid.xs{1};
State{6} = grid.xs{2};
State{7} = tether_diff;

[distance_grid, sol_grid,p_C_W_grid, p_kite_W_grid] = sys.getTargetdistance(State, targetDistanceArgs);
%[gOut, dataOut] = proj(grid, distance_grid, ~[0 0 0 0 1 1], initialState(1:4)');
gOut = grid;
dataOut = distance_grid;
figure(8)
plot3(gOut.xs{1}, gOut.xs{2}, zeros(size(gOut.xs{1})), 'kh'); hold on
h = surf(gOut.xs{1}, gOut.xs{2}, reshape(dataOut, size(gOut.xs{1})));
h.EdgeColor = 'none'; 
h.FaceLighting = 'phong';
