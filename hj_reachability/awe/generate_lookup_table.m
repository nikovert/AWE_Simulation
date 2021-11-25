function [alpha_table, mu_table ,I_table] = generate_lookup_table(g, data, tau, dynSys, extraArgs)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if nargin < 5
    extraArgs = [];
end
   
if isfield(extraArgs, 'newgrid')  && size(data, g.dim+1) == 1
    data = regrid(data, g, extraArgs.newgrid);
    g = extraArgs.newgrid;
end

tauLength = length(tau);
clns = repmat({':'}, 1, g.dim);
derivFunc = @upwindFirstWENO5;

N = g.shape;
grid_min = g.min;
grid_max = g.max;    
dx       = g.dx;
for i = 1:g.dim
  if isfield(g, 'bdry') && isequal(g.bdry{i}, @addGhostPeriodic)
      N(i) = N(i) + 1;
      grid_max(i) = grid_max(i) + dx(i);
  end
  switch i
      case {1, 5, 6}
          dx(i) = dx(i) * dynSys.a0;
          grid_min(i) = grid_min(i) * dynSys.a0;
          grid_max(i) = grid_max(i) * dynSys.a0;
      case 3
          dx(i) = dx(i) * dynSys.h0;
          grid_min(i) = grid_min(i) * dynSys.h0;
          grid_max(i) = grid_max(i) * dynSys.h0;
      case 4
          dx(i) = dx(i) * dynSys.v0;
          grid_min(i) = grid_min(i) * dynSys.v0;
          grid_max(i) = grid_max(i) * dynSys.v0;
  end
          
end

alpha_table = zeros([N tauLength], 'single');
mu_table    = zeros([N tauLength], 'single');
I_table     = zeros([N tauLength], 'uint8');

alpha_options = dynSys.alpha_options;
mu_options = dynSys.mu_options;
values = zeros(length(g.xs{4}(:)), alpha_options*mu_options, 'single');

mu_max = dynSys.mu_max;
alpha_max = dynSys.alpha_max;

for t = 1:tauLength
    BRS_at_t = data(clns{:}, t);
    Deriv = computeGradients(g, BRS_at_t, true(g.dim, 1), derivFunc);
    
    
    q4 = round(Deriv{4}(:), 6);
    q5 = round(Deriv{5}(:), 6)/dynSys.a0;
    q6 = round(Deriv{6}(:), 6)/dynSys.a0;
    va      = g.xs{4}(:);
    gamma_a = g.xs{6}(:) * dynSys.a0;

    fun = @(u) 4.69463 .* q4 .* va .* ((-0.0110857 + 0.183225 .* u(1) +             u(1).^2) .* cos(u(1)) +...
                                   (-0.204926  - 1.98344  .* u(1) + 2.24927  .* u(1).^2) .* sin(u(1))) - ...
                       10.5595 .* ((-0.0911077 - 0.881814 .* u(1) +             u(1).^2) .* cos(u(1)) + ...
                                    (0.0049286 - 0.08146  .* u(1) - 0.444589 .* u(1).^2) .* sin(u(1))) .* ...
                                    (q6 .* cos(u(2)) + q5 .* sec(gamma_a) .* sin(u(2)));

    for m = 1:mu_options
        for a = 1:alpha_options
            alpha = -alpha_max  + (a-1)/(alpha_options-1) * 2 * alpha_max;
            mu    = -mu_max     + (m-1)/(mu_options-1)    * 2 * mu_max;
            values(:, (m-1)*alpha_options + a) = single(fun([alpha, mu]));
        end
    end

    [~, I] = min(values,[],2);
    [I1,I2] = ind2sub([alpha_options, mu_options], I);
    
    alpha = -alpha_max + (I1-1)/(alpha_options-1) * 2 * alpha_max;
    mu    = -mu_max    + (I2-1)/(mu_options-1)    * 2 * mu_max;
    
    alpha = reshape(alpha, size(g.xs{1}));
    mu = reshape(mu, size(g.xs{1}));
    I = reshape(I, size(g.xs{1}));
    
    %uOpt = dynSys.optCtrl(0, g.xs, Deriv);
    
    [~, alpha] = augmentPeriodicData(g, alpha);
    [~, mu] = augmentPeriodicData(g, mu);
    [~, I] = augmentPeriodicData(g, I);
    
    alpha_table(clns{:}, t) = alpha;
    mu_table(clns{:}, t)    = mu;
    I_table(clns{:}, t)     = I;
end

[~, dataf_bool] = augmentPeriodicData(g, int8(data(:,:,:,:,:,:,:,end)<0));

max_tether_diff_dot = dynSys.max_tether_diff_dot;

if isfield(extraArgs, 'saveFile') && extraArgs.saveFile
    %save('tables.mat', 'alpha_table', 'mu_table', 'grid_min', 'grid_max', 'dx', '-v7.3');
    save('tables.mat', 'dataf_bool', 'I_table', 'alpha_options', 'mu_options', 'alpha_max', 'mu_max', 'grid_min', 'grid_max', 'dx', 'max_tether_diff_dot', '-v7.3');
end
end




