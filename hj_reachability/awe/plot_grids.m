clear
close all
%% Grid setup
h0 = 1;
v0 = 1;
a0 = 1;

N        = [       31;   7;       7;      9;        10;         11;     9];
%N        = [       35;   9;       9;      9;        10;         11;     11]; 
grid_min = [     0/a0; -20;  200/h0;  20/v0;    -pi/a0;   -pi/3/a0; -1e-3]; 
grid_max = [  2*pi/a0;  20;  600/h0;  40/v0;     pi/a0;    pi/3/a0;  7e-3];
pdDims   = [1 5];
process  = true;
%% Parameters
h_tau   = 250;
Lem.a       = 120;
Lem.b       = 200;
Lem.phi0    = 1;
extraArgs.Lem = Lem;

Lem.a       = Lem.a /h_tau;
Lem.b       = Lem.b/h_tau;

long = @(s) Lem.b .* sin(s) ./( 1+(Lem.a/Lem.b .* cos(s)).^2 );
lat  = @(s) Lem.a .* sin(s).*cos(s) ./ ( 1+(Lem.a/Lem.b .* cos(s)).^2 );
M_WP = [cos(Lem.phi0),0, -sin(Lem.phi0);0, 1, 0; sin(Lem.phi0),0, cos(Lem.phi0)];

s_curve= 0:0.001:2*pi;
[x,y,z] = sph2cart(long(s_curve),lat(s_curve),h_tau);

points = M_WP * [x;y;z];
%% Figure setup
h1 = figure('WindowState','maximized');
% Create axes
axes1 = axes('Parent',h1);
axis off
hold(axes1,'on');
view(axes1,[108.903466695156 14.0976794597227]);
axis([0, h_tau, -h_tau, h_tau, 0, h_tau])
xlabel('X','Visible','off')
ylabel('Y','Visible','off')
zlabel('Z','Visible','off')

h2 = figure('WindowState','maximized');
% Create axes
axes2 = axes('Parent',h2);
axis off
hold(axes2,'on');
view(axes2,[108.903466695156 14.0976794597227]);
axis([0, h_tau, -h_tau, h_tau, 0, h_tau])
xlabel('X','Visible','off')
ylabel('Y','Visible','off')
zlabel('Z','Visible','off')
%% Plot curve in Fig 1
set(0,'CurrentFigure',h1)

scatter3(points(1,:), points(2,:), points(3,:), 20)
[X,Y,Z] = sphere(55);
surf(X*h_tau,Y*h_tau,Z*h_tau, 'FaceAlpha', 0.15, 'EdgeColor', [0.741176470588235 0.733333333333333 0.733333333333333], ...
    'FaceColor','none','LineStyle', '--', 'Parent',axes1)
line([0 0], [0 0], [h_tau 0],'Color','red','LineStyle','--', 'Parent',axes1,'LineWidth', 2)
line([0 0], [h_tau -h_tau], [0 0],'Color','red','LineStyle','--', 'Parent',axes1,'LineWidth', 2)
line([h_tau 0], [0 0], [0 0],'Color','red','LineStyle','--', 'Parent',axes1,'LineWidth', 2)

%% Plot curve in Fig 2
set(0,'CurrentFigure',h2)

scatter3(points(1,:), points(2,:), points(3,:), 20)
[X,Y,Z] = sphere(55);
surf(X*h_tau,Y*h_tau,Z*h_tau, 'FaceAlpha', 0.15, 'EdgeColor', [0.741176470588235 0.733333333333333 0.733333333333333], ...
    'FaceColor','none','LineStyle', '--', 'Parent',axes2)
line([0 0], [0 0], [h_tau 0],'Color','red','LineStyle','--', 'Parent',axes2,'LineWidth', 2)
line([0 0], [h_tau -h_tau], [0 0],'Color','red','LineStyle','--', 'Parent',axes2,'LineWidth', 2)
line([h_tau 0], [0 0], [0 0],'Color','red','LineStyle','--', 'Parent',axes2,'LineWidth', 2)

%% Plot sol/sigma grid
set(0,'CurrentFigure',h1)
s = linspace(grid_min(1), grid_max(1), N(1));
sigma = linspace(grid_min(2), grid_max(2), N(2));
[S, SIGMA] = meshgrid(s,sigma);

[long, lat] = getLongLat(S, SIGMA, h_tau, extraArgs);

[x,y,z] = sph2cart(long(:),lat(:),h_tau);
scatter3(x, y, z, 40, 'ko', 'filled')
%% Plot long/lat grid
set(0,'CurrentFigure',h2)
long = linspace(min(long(:)), max(long(:)), N(1));
lat = linspace(min(lat(:)), max(lat(:)), N(2));
[LONG, LAT] = meshgrid(long,lat);

[x,y,z] = sph2cart(LONG(:),LAT(:),h_tau);
scatter3(x, y, z, 40,'ko' , 'filled')