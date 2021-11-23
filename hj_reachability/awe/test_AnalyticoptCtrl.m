%% Add relevant files to path
addpath('../')
addToolbox;
    
load('optCtrl_workspace.mat');
va = va(1);
gamma_a = gamma_a(1);
q4 = q4(1);
q5 = q5(1);
q6 = q6(1);
%% Define alpha part
a0 = -0.05204341693058464;
a1 = 0.8601751285622818;
a2 = 4.694628750000001;

b0 = -0.9620506210503925;
b1 = -9.31149481764243;
b2 = 10.55948284375565;

a = @(alpha) a0 + a1 * alpha + a2 * alpha.^2;
dadt = @(alpha) a1 + 2 * a2 * alpha;

b = @(alpha) b0 + b1 * alpha + b2 * alpha.^2;
dbdt = @(alpha) b1 + 2 * b2 * alpha; 

f1 = @(alpha) a(alpha) .* cos(alpha) + b(alpha) .* sin(alpha);
df1dt = @(alpha) dadt(alpha) .* cos(alpha) + dbdt(alpha) .* sin(alpha) - ...
        a(alpha) .* sin(alpha) + b(alpha) .* cos(alpha);

f2 = @(alpha) b(alpha) .* cos(alpha) - a(alpha) .* sin(alpha);
df2dt = @(alpha) dbdt(alpha) .* cos(alpha) - dadt(alpha) .* sin(alpha) - ...
        b(alpha) .* sin(alpha) - a(alpha) .* cos(alpha);
    
h1 = figure(1);
grid on; hold on;
plot(-pi:0.1:pi, f1(-pi:0.1:pi))
plot(-pi:0.1:pi, f2(-pi:0.1:pi))

%% Define mu part
g1 = va .* q4;
g2 = @(mu)(q6.^2 + q5.^2 .* sec(gamma_a).^2) .* cos(mu - atan2(q5,(cos(gamma_a) .* q6)));

h2 = figure(2);
grid on; hold on;
plot(-pi:0.1:pi, g2(-pi:0.1:pi))
%% See if we can find the minimum
H = @(alpha, mu) g1.* f1(alpha) + g2(mu) .* f2(alpha);

% Case 1
mu = atan2(q5,(cos(gamma_a) .* q6));
alpha1 = fminbnd(@(alpha) g1.* f1(alpha) + g2(mu) .* f2(alpha),-obj.alpha_max,obj.alpha_max);
H(alpha1, mu)

% Case 2
mu = atan2(q5,(cos(gamma_a) .* q6)) + pi/2;
alpha2 = fminbnd(@(alpha) g1.* f1(alpha) + g2(mu) .* f2(alpha),-obj.alpha_max,obj.alpha_max);
H(alpha2, mu)

% Case 3
mu = atan2(q5,(cos(gamma_a) .* q6)) + pi;
alpha3 = fminbnd(@(alpha) g1.* f1(alpha) + g2(mu) .* f2(alpha),-obj.alpha_max,obj.alpha_max);
H(alpha3, mu)

%%
h3 = figure(3);
grid on; hold on;
[ALPHA, MU] = meshgrid(-obj.alpha_max:0.1:obj.alpha_max, -obj.mu_max:0.1:obj.mu_max);
surf(ALPHA, MU, F(ALPHA, MU))
view(45,45)