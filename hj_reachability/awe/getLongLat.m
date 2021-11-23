function [long, lat, t_tau, t_rot_tau, t_W, t_rot_W] = getLongLat(s, sigma, h_tau, extraArgs)

if nargin < 4
    extraArgs = [];
end

if ~isfield(extraArgs, 'Lem')
    Lem.a = 120;
    Lem.b = 200;
    Lem.phi0 = 1;
else
    Lem = extraArgs.Lem;
end

if ~isfield(extraArgs, 'visualize_steps')
    visualize_steps = false;
else
    visualize_steps = extraArgs.visualize_steps;
end

if length(s) > 1 && length(sigma) > 1
    assert(length(s) == length(sigma))
end

if ~isfield(extraArgs, 'direction')
    direction = 1;
else
    direction = extraArgs.direction;
end

M_WP = {cos(Lem.phi0),0, -sin(Lem.phi0);0, 1, 0; sin(Lem.phi0),0, cos(Lem.phi0)};
Lem.a = Lem.a ./ h_tau;
Lem.b = Lem.b ./ h_tau;

%% calculate long_s and lat_s
a = Lem.a;
b = Lem.b;
L = { b .* sin(s) ./( 1+(a./b .* cos(s)).^2 );
          a .* sin(s).*cos(s) ./ ( 1+(a./b .* cos(s)).^2 ) } ;

dLds = { ( b.^3 .* cos(s).*(2*a.^2-a.^2 .* cos(s).^2+b.^2)./(a.^2 .* cos(s).^2+b.^2).^2 );
( (cos(s).^2 .* (a.^3 .* b.^2+2*a .* b.^4) - a .* b.^4)./(a.^2 .* cos(s).^2+b.^2).^2 )};  

s_lambda = sin( L{1}  ); 
s_phi = sin( L{2} );
c_lambda = cos( L{1} ); 
c_phi = cos( L{2} ); 

dqdlambda = {-s_lambda.*c_phi; c_lambda.*c_phi; zeros(size(s_phi))};
dqdphi = {-c_lambda.*s_phi; -s_lambda.*s_phi; c_phi};

t_P = {direction*(dqdlambda{1} .* dLds{1} + dqdphi{1} .* dLds{2}); ...
         direction*(dqdlambda{2} .* dLds{1} + dqdphi{2} .* dLds{2}); ...
         direction*(dqdlambda{3} .* dLds{1} + dqdphi{3} .* dLds{2})};      
      

t_W = mult_cellMatrix(M_WP,t_P); 

[pos_C_P_x,pos_C_P_y,pos_C_P_z] = sph2cart(L{1},L{2},h_tau);
p_C_P = {pos_C_P_x;pos_C_P_y;pos_C_P_z};
    
p_C_W = mult_cellMatrix(M_WP , p_C_P);
[long_C_W, lat_C_W, ~] = cart2sph(p_C_W{1}, p_C_W{2}, p_C_W{3});

%% move out by delta
M_tauW = {-sin(lat_C_W).*cos(long_C_W), -sin(lat_C_W).*sin(long_C_W),  cos(lat_C_W);
          -sin(long_C_W)          , cos(long_C_W)           ,         0;
          -cos(lat_C_W).*cos(long_C_W), -cos(lat_C_W).*sin(long_C_W), -sin(lat_C_W)};
      
M_Wtau = transpose(M_tauW);

t_tau = mult_cellMatrix(M_tauW , t_W);
M_rotation = {cos(direction*pi/2), -sin(direction*pi/2),  0; ...
              sin(direction*pi/2), cos(direction*pi/2), 0; ...
              0, 0, 1};
      
t_rot_tau = mult_cellMatrix(M_rotation , t_tau);
t_rot_W = mult_cellMatrix(M_Wtau , t_rot_tau);

[long_dot, lat_dot, ~] = vel_cart2sph(p_C_W{1}, p_C_W{2}, p_C_W{3}, t_rot_W{1}./norm_cellVec(t_rot_W),t_rot_W{2}./norm_cellVec(t_rot_W),t_rot_W{3}./norm_cellVec(t_rot_W));
X_initial = {long_C_W; long_dot; lat_C_W; lat_dot};
%[~,X]=ode45(@cir,[0,sigma],X_initial);

N = 10;
h = sigma/N;
X = X_initial;
for k=1:N
    X = add_cellMatrix(X , mult_cellMatrix({h}, cir(0,X)));
end
u=X{1};
v=X{3};

[x_geo, y_geo, z_geo] = sph2cart(u, v, h_tau);

p_kite_W = {x_geo, y_geo, z_geo};

long = atan(p_kite_W{2}./p_kite_W{1});
lat = asin(p_kite_W{3}./h_tau);
% if abs(norm_cellVec(p_kite_W) - h_tau) > 0.001
%     warning('h_tau does not match')
% end
%[p_kite_W{1}, p_kite_W{2}, p_kite_W{3}] = sph2cart(long, lat, h_tau);

% t_tau_normalised = element_div_cellMatrix(t_tau,(norm_cellVec(t_tau)));
% % assert(abs(sin(acos(t_tau_normalised{2})) - t_tau_normalised{1}) < 1e-8)
% % assert(abs(cos(asin(t_tau_normalised{1})) - t_tau_normalised{2}) < 1e-8)
% theta = asin(t_tau_normalised{1});
% M_Ltau = {cos(theta), -sin(theta),  0; ...
%           sin(theta), cos(theta), 0; ...
%           0, 0, 1};
%       warning('check rotation direction of M_Ltau')
          
if visualize_steps
    sol_spread = 0:0.01:2*pi;
    i = length(h_tau);
    long_path = Lem.b(i) .* sin(sol_spread) ./ ( 1+(Lem.a(i) ./ Lem.b(i) .* cos(sol_spread)).^2 );
    lat_path  = Lem.a(i) .* sin(sol_spread) .* cos(sol_spread) ./ ( 1+(Lem.a(i) ./Lem.b(i) .* cos(sol_spread)).^2 ) ;
    [X,Y,Z] = sph2cart(long_path,lat_path,h_tau(1));
    points = cell2mat(M_WP) * [X;Y;Z];
    scatter3(points(1,:),points(2,:),points(3,:),20)
    hold on
    
    scale_factor = 100;
    % visualize velocity at tangent of desired point
    quiver3(p_C_W{1},p_C_W{2},p_C_W{3},scale_factor*t_W{1},scale_factor*t_W{2},scale_factor*t_W{3})
    quiver3(p_C_W{1},p_C_W{2},p_C_W{3},translation_vec{1},translation_vec{2},translation_vec{3})
    plot3(p_kite_W{1}, p_kite_W{2}, p_kite_W{3}, 'kh')
end

end

function xp = cir(t,x)
    xp = cell(4,1);
    xp{1} = x{2};
    xp{2} = 2*tan(x{3}).*x{2}.*x{4};
    xp{3} = x{4};
    xp{4} = -sin(x{3}).*cos(x{3}).*x{2}.^2;
end