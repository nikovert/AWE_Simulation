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

function [alpha_table, mu_table ,I_table] = generate_lookup_table(g, data, tau, dynSys, extraArgs)
%GENERATE_LOOKUP_TABLE generates the looup table used during Simulation
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
d1_table    = zeros([N tauLength], 'single');
d2_table    = zeros([N tauLength], 'single');
d3_table    = zeros([N tauLength], 'single');
d4_table    = zeros([N tauLength], 'single');
I_table     = zeros([N tauLength], 'uint8');

alpha_options = dynSys.alpha_options;
mu_options = dynSys.mu_options;
values = zeros(length(g.xs{4}(:)), alpha_options*mu_options, 'single');

mu_max = dynSys.mu_max;
mu_min = dynSys.mu_min;
alpha_max = dynSys.alpha_max;
alpha_min = dynSys.alpha_min;

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
            alpha = alpha_min  + (a-1)/(alpha_options-1) * (alpha_max-alpha_min);
            mu    = mu_min     + (m-1)/(mu_options-1)    * (mu_max-mu_min);
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

    % disturbances
    dOpt = dynSys.optDstb(t, g.xs, Deriv);
    [~, dOpt{1}] = augmentPeriodicData(g, dOpt{1});
    [~, dOpt{2}] = augmentPeriodicData(g, dOpt{2});
    [~, dOpt{3}] = augmentPeriodicData(g, dOpt{3});
    [~, dOpt{4}] = augmentPeriodicData(g, dOpt{4});
    d1_table(clns{:}, t)    = single(dOpt{1});
    d2_table(clns{:}, t)    = single(dOpt{2});
    d3_table(clns{:}, t)    = single(dOpt{3});
    d4_table(clns{:}, t)    = single(dOpt{4});
end

[~, dataf_bool] = augmentPeriodicData(g, int8(data(:,:,:,:,:,:,:,end)<0));

max_tether_diff_dot = dynSys.max_tether_diff_dot;

if isfield(extraArgs, 'saveFile') && extraArgs.saveFile
    save('tables.mat', 'dataf_bool', 'I_table', 'alpha_options', ...
        'mu_options', 'alpha_max', 'alpha_min', 'mu_max', 'mu_min', ...
        'grid_min', 'grid_max', 'dx', 'max_tether_diff_dot',...
        'd1_table', 'd2_table', 'd3_table', 'd4_table', '-v7.3');
end
end




