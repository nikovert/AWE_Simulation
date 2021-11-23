function [long, lat, t_tau, t_rot_tau, t_W, t_rot_W] = getLongLat(s, sigma, h_tau, PHI_BOOTH, params)

Lem.a = params.a_booth;
Lem.b = params.b_booth;
Lem.phi0 = PHI_BOOTH;

direction = params.direction;

M_WP = [cos(Lem.phi0), 0, -sin(Lem.phi0);...
                    0, 1,              0;...
        sin(Lem.phi0), 0,  cos(Lem.phi0)];
Lem.a = Lem.a ./ h_tau;
Lem.b = Lem.b ./ h_tau;

%% calculate long_s and lat_s
[t_P,~,L, ~] = getBoothInfos2(s,Lem, direction);
t_W = M_WP*t_P; 

[pos_C_P_x,pos_C_P_y,pos_C_P_z] = sph2cart(L(1),L(2),h_tau);
p_C_P = [pos_C_P_x;pos_C_P_y;pos_C_P_z];
    
p_C_W = M_WP* p_C_P;
[long_C_W, lat_C_W, ~] = cart2sph(p_C_W(1), p_C_W(2), p_C_W(3));

%% move out by delta
M_tauW = [-sin(lat_C_W).*cos(long_C_W), -sin(lat_C_W).*sin(long_C_W),  cos(lat_C_W);
          -sin(long_C_W)          , cos(long_C_W)           ,         0;
          -cos(lat_C_W).*cos(long_C_W), -cos(lat_C_W).*sin(long_C_W), -sin(lat_C_W)];
      
M_Wtau = M_tauW';

t_tau = M_tauW * t_W;
M_rotation = [cos(direction * pi/2), -sin(direction * pi/2),  0; ...
              sin(direction * pi/2), cos(direction * pi/2), 0; ...
              0, 0, 1];
      
t_rot_tau = M_rotation * t_tau;
t_rot_W = M_Wtau * t_rot_tau;

[long_dot, lat_dot, ~] = vel_cart2sph(p_C_W(1), p_C_W(2), p_C_W(3), t_rot_W(1)/norm(t_rot_W),t_rot_W(2)/norm(t_rot_W),t_rot_W(3)/norm(t_rot_W));
X_initial = [long_C_W, long_dot, lat_C_W, lat_dot];
%[~,X]=ode45(@cir,[0,sigma],X_initial);

N = 100;
h = sigma/N;
X = zeros(length(X_initial),N+1);
X(:,1) = X_initial;
for k=1:N
    X(:,k+1) = X(:,k) + h * cir(0,X(:,k));
end

u=X(1,end);
v=X(3,end);

[x_geo, y_geo, z_geo] = sph2cart(u, v, h_tau);

p_kite_W = [x_geo, y_geo, z_geo];

if ~any(isreal(p_kite_W)) || ~isreal(h_tau)
    1;
elseif (p_kite_W(3)/h_tau) > 1
    1;
elseif (p_kite_W(3)/h_tau) < -1
    1;
end
long = atan(p_kite_W(2)./p_kite_W(1));
lat = asin(max(-1,min(1,p_kite_W(3)./h_tau)));
end

function xp = cir(t,x)
    xp = zeros(4,1);
    xp(1) = x(2);
    xp(2) = 2*tan(x(3))*x(2)*x(4);
    xp(3) = x(4);
    xp(4) = -sin(x(3))*cos(x(3))*x(2)^2;
end