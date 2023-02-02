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

function uOpt = optCtrl(obj, t, x, deriv, uMode, ~)
% determine the optimal input (i.e. that maximises/minimises the Hamiltonian)

if nargin < 5
  uMode = 'min';
end

if ~(strcmp(uMode, 'max') || strcmp(uMode, 'min'))
  error('uMode must be ''max'' or ''min''!')
end

num = false;
if ~iscell(x)
  x = num2cell(x);
  num = true;
end

if ~iscell(deriv)
  deriv = num2cell(deriv);
end

q4 = round(deriv{4}(:), 6);
q5 = round(deriv{5}(:), 6)/obj.a0;
q6 = round(deriv{6}(:), 6)/obj.a0;

va      = x{4}(:);
gamma_a = x{6}(:) * obj.a0;
% if ~isempty(obj.F_rest) && ~num
%     F_rest = obj.F_rest;
% elseif num
%    [~, F_rest] = obj.dynamics(t, cell2mat(x), {0;0});
% else
%     [~, F_rest] = obj.dynamics(t, x, {0;0});
% end
% Frest = {F_rest{1}(:); F_rest{2}(:); F_rest{3}(:)};

% fun_old = @(u) ((4.69463 .* (q4./obj.v0) .* (va*obj.v0) .* (0.213009 .* Frest{1}+(-0.0110857+0.183225 .* u(1)+u(1).^2) .* (va*obj.v0).^2 .* cos(u(1)) ...
%             +(-0.204926-1.98344 .* u(1)+2.24927 .* u(1).^2) .* (va*obj.v0).^2 .* sin(u(1)))+Frest{2} .* (q5 .* cos(u(2)) .* sec(gamma_a) ...
%             -q6 .* sin(u(2)))- (Frest{3}+(-0.962051-9.31149 .* u(1)+10.5595 .* u(1).^2) .* (va*obj.v0).^2 .* cos(u(1))+(0.0520434 ...
%             -0.860175 .* u(1)-4.69463 .* u(1).^2) .* (va*obj.v0).^2 .* sin(u(1))) .* (q6 .* cos(u(2))+q5 .* sec(gamma_a) .* sin(u(2))))./(va*obj.v0));

% fun = @(u) 4.69463 .* q4 .* va .* ((-0.0110857 + 0.183225 .* u(1) +             u(1).^2) .* cos(u(1)) +...
%                                    (-0.204926  - 1.98344  .* u(1) + 2.24927  .* u(1).^2) .* sin(u(1))) - ...
%                        10.5595 .* ((-0.0911077 - 0.881814 .* u(1) +             u(1).^2) .* cos(u(1)) + ...
%                                     (0.0049286 - 0.08146  .* u(1) - 0.444589 .* u(1).^2) .* sin(u(1))) .* ...
%                                     (q6 .* cos(u(2)) + q5 .* sec(gamma_a) .* sin(u(2)));
    
alpha = opt_alpha(va, obj);
alpha = alpha(:);
a = 4.69463 .* q4 .* va .* ((-0.0110857 + 0.183225 .* alpha +             alpha.^2) .* cos(alpha) +...
                                   (-0.204926  - 1.98344  .* alpha + 2.24927  .* alpha.^2) .* sin(alpha)) - ...
                       10.5595 .* ((-0.0911077 - 0.881814 .* alpha +             alpha.^2) .* cos(alpha) + ...
                                    (0.0049286 - 0.08146  .* alpha - 0.444589 .* alpha.^2) .* sin(alpha)) .* q6;
b = 4.69463 .* q4 .* va .* ((-0.0110857 + 0.183225 .* alpha +             alpha.^2) .* cos(alpha) +...
                                   (-0.204926  - 1.98344  .* alpha + 2.24927  .* alpha.^2) .* sin(alpha)) - ...
                       10.5595 .* ((-0.0911077 - 0.881814 .* alpha +             alpha.^2) .* cos(alpha) + ...
                                    (0.0049286 - 0.08146  .* alpha - 0.444589 .* alpha.^2) .* sin(alpha)) .* q5 .* sec(gamma_a);

mu = atan2(b,a);
if strcmp(uMode, 'min')
    mu = mu + pi;
end

% fun = @(mu) 4.69463 .* q4 .* va .* ((-0.0110857 + 0.183225 .* alpha +             alpha.^2) .* cos(alpha) +...
%                                    (-0.204926  - 1.98344  .* alpha + 2.24927  .* alpha.^2) .* sin(alpha)) - ...
%                        10.5595 .* ((-0.0911077 - 0.881814 .* alpha +             alpha.^2) .* cos(alpha) + ...
%                                     (0.0049286 - 0.08146  .* alpha - 0.444589 .* alpha.^2) .* sin(alpha)) .* ...
%                                     (q6 .* cos(mu) + q5 .* sec(gamma_a) .* sin(mu));
% 
% alpha_options = obj.alpha_options;
% mu_options = obj.mu_options;
% values = zeros(length(single(fun([0, 0]))), alpha_options*mu_options, 'single');
% for m = 1:mu_options
%     for a = 1:alpha_options
%         alpha = obj.alpha_min + (a-1)/(alpha_options-1) * 2 * obj.alpha_max;
%         mu    = obj.mu_min + (m-1)/(mu_options-1) * 2 * obj.mu_max;
%         values(:,(m-1)*alpha_options + a) = single(fun([alpha, mu]));
%     end
% end
% % values_old = [fun_old([obj.alpha_min,obj.mu_min])    fun_old([obj.alpha_min/2,obj.mu_min])    fun_old([0,obj.mu_min])    fun_old([obj.alpha_max/2,obj.mu_min])  fun_old([obj.alpha_max,obj.mu_min]) ...
% %               fun_old([obj.alpha_min,obj.mu_min/2])    fun_old([obj.alpha_min/2,obj.mu_min/2])    fun_old([0,obj.mu_min/2])    fun_old([obj.alpha_max/2,obj.mu_min/2])  fun_old([obj.alpha_max,obj.mu_min/2]) ...
% %               fun_old([obj.alpha_min,0])              fun_old([obj.alpha_min/2,0])              fun_old([0,0])              fun_old([obj.alpha_max/2,0])            fun_old([obj.alpha_max,0]) ...
% %               fun_old([obj.alpha_min,obj.mu_max/2])    fun_old([obj.alpha_min/2,obj.mu_max/2])    fun_old([0,obj.mu_max/2])    fun_old([obj.alpha_max/2,obj.mu_max/2])  fun_old([obj.alpha_max,obj.mu_max/2]) ...
% %               fun_old([obj.alpha_min,obj.mu_max])     fun_old([obj.alpha_min/2,obj.mu_max])     fun_old([0,obj.mu_max])     fun_old([obj.alpha_max/2,obj.mu_max])   fun_old([obj.alpha_max,obj.mu_max])];       
% 
% % assert(mu_options == 9 && alpha_options == 5);
% % values = [fun([obj.alpha_min,    obj.mu_min])    fun([obj.alpha_min/2,    obj.mu_min])    fun([0,    obj.mu_min])    fun([obj.alpha_max/2,    obj.mu_min])    fun([obj.alpha_max,    obj.mu_min])      ...
% %           fun([obj.alpha_min,-3*obj.mu_max/4])    fun([obj.alpha_min/2,-3*obj.mu_max/4])    fun([0,-3*obj.mu_max/4])    fun([obj.alpha_max/2,-3*obj.mu_max/4])    fun([obj.alpha_max,-3*obj.mu_max/4])      ...
% %           fun([obj.alpha_min,  obj.mu_min/2])    fun([obj.alpha_min/2,  obj.mu_min/2])    fun([0,  obj.mu_min/2])    fun([obj.alpha_max/2,  obj.mu_min/2])    fun([obj.alpha_max,  obj.mu_min/2])      ...
% %           fun([obj.alpha_min,  obj.mu_min/4])    fun([obj.alpha_min/2,  obj.mu_min/4])    fun([0,  obj.mu_min/4])    fun([obj.alpha_max/2,  obj.mu_min/4])    fun([obj.alpha_max,  obj.mu_min/4])      ...
% %           fun([obj.alpha_min,              0])    fun([obj.alpha_min/2,              0])    fun([0,              0])    fun([obj.alpha_max/2,              0])    fun([obj.alpha_max,              0])      ...
% %           fun([obj.alpha_min,   obj.mu_max/4])    fun([obj.alpha_min/2,   obj.mu_max/4])    fun([0,   obj.mu_max/4])    fun([obj.alpha_max/2,   obj.mu_max/4])    fun([obj.alpha_max,   obj.mu_max/4])      ...
% %           fun([obj.alpha_min,   obj.mu_max/2])    fun([obj.alpha_min/2,   obj.mu_max/2])    fun([0,   obj.mu_max/2])    fun([obj.alpha_max/2,   obj.mu_max/2])    fun([obj.alpha_max,   obj.mu_max/2])      ...
% %           fun([obj.alpha_min, 3*obj.mu_max/4])    fun([obj.alpha_min/2, 3*obj.mu_max/4])    fun([0, 3*obj.mu_max/4])    fun([obj.alpha_max/2, 3*obj.mu_max/4])    fun([obj.alpha_max, 3*obj.mu_max/4])    ...
% %           fun([obj.alpha_min,     obj.mu_max])    fun([obj.alpha_min/2,     obj.mu_max])    fun([0,     obj.mu_max])    fun([obj.alpha_max/2,     obj.mu_max])    fun([obj.alpha_max,     obj.mu_max])];       
% 
% if strcmp(uMode, 'min')
%     %[~, I] = min(values_old,[],2);
%     [~, I] = min(values,[],2);
% else
%     [~, I] = max(values,[],2);
% end
% [a,m] = ind2sub([alpha_options, mu_options],I);
% alpha = obj.alpha_min + (a-1)/(alpha_options-1) * 2 * obj.alpha_max;
% mu    = obj.mu_min + (m-1)/(mu_options-1) * 2 * obj.mu_max;
alpha = reshape(alpha, size(x{1}));
mu = reshape(mu, size(x{1}));

uOpt = {alpha; mu};
if num
    uOpt = [alpha; mu];
end

end

% Calculate alpha for desired lift force
function alpha = opt_alpha(Va, sys)
    L_rq = sys.Ft_set;
    CL_rq = 2*L_rq / (sys.AIRCRAFT.rho_air*sys.AIRCRAFT.S_wing*Va.^2 );
    alpha = (CL_rq - sys.AIRCRAFT.CL0)/sys.AIRCRAFT.CL_alpha; 
end