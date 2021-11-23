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

fun = @(u) 4.69463 .* q4 .* va .* ((-0.0110857 + 0.183225 .* u(1) +             u(1).^2) .* cos(u(1)) +...
                                   (-0.204926  - 1.98344  .* u(1) + 2.24927  .* u(1).^2) .* sin(u(1))) - ...
                       10.5595 .* ((-0.0911077 - 0.881814 .* u(1) +             u(1).^2) .* cos(u(1)) + ...
                                    (0.0049286 - 0.08146  .* u(1) - 0.444589 .* u(1).^2) .* sin(u(1))) .* ...
                                    (q6 .* cos(u(2)) + q5 .* sec(gamma_a) .* sin(u(2)));
    
alpha_options = obj.alpha_options;
mu_options = obj.mu_options;
assert(mu_options == 9 && alpha_options == 5);
% values = zeros(length(q4), alpha_options*mu_options);
% for m = 1:mu_options
%     for a = 1:alpha_options
%         alpha = -obj.alpha_max + (a-1)/(alpha_options-1) * 2 * obj.alpha_max;
%         mu    = -obj.mu_max + (m-1)/(mu_options-1) * 2 * obj.mu_max;
%         values((m-1)*alpha_options + a) = fun([alpha, mu]);
%     end
% end
% values_old = [fun_old([-obj.alpha_max,-obj.mu_max])    fun_old([-obj.alpha_max/2,-obj.mu_max])    fun_old([0,-obj.mu_max])    fun_old([obj.alpha_max/2,-obj.mu_max])  fun_old([obj.alpha_max,-obj.mu_max]) ...
%               fun_old([-obj.alpha_max,-obj.mu_max/2])    fun_old([-obj.alpha_max/2,-obj.mu_max/2])    fun_old([0,-obj.mu_max/2])    fun_old([obj.alpha_max/2,-obj.mu_max/2])  fun_old([obj.alpha_max,-obj.mu_max/2]) ...
%               fun_old([-obj.alpha_max,0])              fun_old([-obj.alpha_max/2,0])              fun_old([0,0])              fun_old([obj.alpha_max/2,0])            fun_old([obj.alpha_max,0]) ...
%               fun_old([-obj.alpha_max,obj.mu_max/2])    fun_old([-obj.alpha_max/2,obj.mu_max/2])    fun_old([0,obj.mu_max/2])    fun_old([obj.alpha_max/2,obj.mu_max/2])  fun_old([obj.alpha_max,obj.mu_max/2]) ...
%               fun_old([-obj.alpha_max,obj.mu_max])     fun_old([-obj.alpha_max/2,obj.mu_max])     fun_old([0,obj.mu_max])     fun_old([obj.alpha_max/2,obj.mu_max])   fun_old([obj.alpha_max,obj.mu_max])];       

values = [fun([-obj.alpha_max,    -obj.mu_max])    fun([-obj.alpha_max/2,    -obj.mu_max])    fun([0,    -obj.mu_max])    fun([obj.alpha_max/2,    -obj.mu_max])    fun([obj.alpha_max,    -obj.mu_max])      ...
          fun([-obj.alpha_max,-3*obj.mu_max/4])    fun([-obj.alpha_max/2,-3*obj.mu_max/4])    fun([0,-3*obj.mu_max/4])    fun([obj.alpha_max/2,-3*obj.mu_max/4])    fun([obj.alpha_max,-3*obj.mu_max/4])      ...
          fun([-obj.alpha_max,  -obj.mu_max/2])    fun([-obj.alpha_max/2,  -obj.mu_max/2])    fun([0,  -obj.mu_max/2])    fun([obj.alpha_max/2,  -obj.mu_max/2])    fun([obj.alpha_max,  -obj.mu_max/2])      ...
          fun([-obj.alpha_max,  -obj.mu_max/4])    fun([-obj.alpha_max/2,  -obj.mu_max/4])    fun([0,  -obj.mu_max/4])    fun([obj.alpha_max/2,  -obj.mu_max/4])    fun([obj.alpha_max,  -obj.mu_max/4])      ...
          fun([-obj.alpha_max,              0])    fun([-obj.alpha_max/2,              0])    fun([0,              0])    fun([obj.alpha_max/2,              0])    fun([obj.alpha_max,              0])      ...
          fun([-obj.alpha_max,   obj.mu_max/4])    fun([-obj.alpha_max/2,   obj.mu_max/4])    fun([0,   obj.mu_max/4])    fun([obj.alpha_max/2,   obj.mu_max/4])    fun([obj.alpha_max,   obj.mu_max/4])      ...
          fun([-obj.alpha_max,   obj.mu_max/2])    fun([-obj.alpha_max/2,   obj.mu_max/2])    fun([0,   obj.mu_max/2])    fun([obj.alpha_max/2,   obj.mu_max/2])    fun([obj.alpha_max,   obj.mu_max/2])      ...
          fun([-obj.alpha_max, 3*obj.mu_max/4])    fun([-obj.alpha_max/2, 3*obj.mu_max/4])    fun([0, 3*obj.mu_max/4])    fun([obj.alpha_max/2, 3*obj.mu_max/4])    fun([obj.alpha_max, 3*obj.mu_max/4])    ...
          fun([-obj.alpha_max,     obj.mu_max])    fun([-obj.alpha_max/2,     obj.mu_max])    fun([0,     obj.mu_max])    fun([obj.alpha_max/2,     obj.mu_max])    fun([obj.alpha_max,     obj.mu_max])];       

if strcmp(uMode, 'min')
    %[~, I] = min(values_old,[],2);
    [~, I] = min(values,[],2);
else
    [~, I] = max(values,[],2);
end
[a,m] = ind2sub([alpha_options, mu_options],I);
alpha = -obj.alpha_max + (a-1)/(alpha_options-1) * 2 * obj.alpha_max;
mu    = -obj.mu_max + (m-1)/(mu_options-1) * 2 * obj.mu_max;
alpha = reshape(alpha, size(x{1}));
mu = reshape(mu, size(x{1}));

uOpt = {alpha; mu};
if num
    uOpt = [alpha; mu];
end

end