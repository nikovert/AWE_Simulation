function [f_kite_G, final_seg_diff_dot, max_F_tether,ground_tether_force] = tether_forces(obj, x_W)
%TETHER_FORCES Calculates the tether forces for a given state x
% Variables of interesst:
% p_vec: Tether particle position (wind frame)
% v_vec: Tether particle velocities (wind frame)
% c0: intial spring constant of a tether segment
% d0: initial damping constant of a tether segment
% rho_t: tether mass density (i.e. mass per m)
% d_tether: tether diameter
% CD_tether: tether drag coefficient
% l_s: tether segment length
% g_vec: gravity vector
%
% obj.ENVMT     - environmental variable
% obj.AIRCRAFT  - Aircarft variables
% T         - tether properties 
%
%% Define states and inputs (remove this step later for speed up)
long        = x_W{1};
lat         = x_W{2};
h_tau       = x_W{3};
long_dot    = x_W{4};
lat_dot     = x_W{5};
h_tau_dot   = x_W{6};
tether_diff = x_W{7};

c0 = obj.T.c0;
d0 = obj.T.d0;
rho_t = obj.AIRCRAFT.rho_tether;% T.rho_t; 
d_tether = obj.AIRCRAFT.d_tether;%  T.d_tether;
CD_tether = obj.AIRCRAFT.Cd_tether;% T.CD_tether;

n_t_p = obj.T.np; % amount of particles

% Adapt stiffness and damping with respect to current segment length
c = c0/60; 
d = d0/60;

% Mass of tether segments with  respect to current segment length
m = h_tau/(n_t_p + 1) * rho_t;

p_vec = cell(n_t_p*3, 1);
v_vec = cell(n_t_p*3, 1);
for i = 1:3:n_t_p*3
   [pos_x,pos_y,pos_z] = sph2cart(long,lat,ceil(i/3) * h_tau/(n_t_p + 1));
   [vel_x,vel_y,vel_z] = vel_sph2cart(long,lat,ceil(i/3) * h_tau/(n_t_p + 1), long_dot,lat_dot,h_tau_dot);
   p_vec{i}   = pos_x;
   p_vec{i+1} = pos_y;
   p_vec{i+2} = pos_z;
   v_vec{i}   = vel_x;
   v_vec{i+1} = vel_y;
   v_vec{i+2} = vel_z;
end
pos_G_W   = cell(3,1);
vel_k_G_W = cell(3,1);

[pos_G_W_x,pos_G_W_y,pos_G_W_z] = sph2cart(long,lat,h_tau);
pos_G_W{1} = pos_G_W_x;
pos_G_W{2} = pos_G_W_y;
pos_G_W{3} = pos_G_W_z;

[vel_k_G_W_x,vel_k_G_W_y,vel_k_G_W_z] = vel_sph2cart(long,lat,h_tau, long_dot,lat_dot,h_tau_dot);
vel_k_G_W{1} = vel_k_G_W_x;
vel_k_G_W{2} = vel_k_G_W_y;
vel_k_G_W{3} = vel_k_G_W_z;

p_norm = cell(n_t_p,1);

% Vector with distances between particles
p_d = sub_cellMatrix(p_vec,[0;0;0; p_vec(1:3*n_t_p-3)]);

% Vector with velocity differences
v_d = sub_cellMatrix(v_vec,[0;0;0; v_vec(1:3*n_t_p-3)]);

% calculate absolute distances between points
for iter = 1:3:n_t_p*3
    p_norm{ceil(iter/3)} = sqrt(p_d{iter}.^2 + p_d{iter+1}.^2 + p_d{iter+2}.^2);
end

% Gravity vector in W frame
Fg= {0;0;m*obj.ENVMT.g};

% vel_W_W_mat = wind velocity for the differnt particles (3 x nr_segments)
z_0 = 0.15;
v_w_u = cell(n_t_p, 1);
for i = 1:n_t_p
    v_w_u{i} = obj.base_windspeed * log(3.281 * p_vec{3} * i/z_0)/log(20/z_0);
end

%% Tether dynamics
% Initialize matrices
fsd_tot_mat = cell(3,n_t_p);
kite_tether_force = cell(3,n_t_p);
D_s_mat = cell(3,n_t_p);

tmp_p = sub_cellMatrix(pos_G_W,p_vec(3*(n_t_p) - 2 : 3*(n_t_p)));
f_kite_G = obj.calcSpringDamperForce( c0, norm_cellVec(tmp_p), tether_diff, d0,tmp_p ,...
    sub_cellMatrix(vel_k_G_W,v_vec(3*(n_t_p) - 2 : 3*(n_t_p))));

final_seg_diff_dot = element_div_cellMatrix(mult_cellMatrix(transpose_cellMatrix(tmp_p), ...
    sub_cellMatrix(vel_k_G_W,v_vec(3*(n_t_p) - 2 : 3*(n_t_p)))) , norm_cellVec(tmp_p));

if nargout < 3
    return
end
    
if false
    X = pos_G_W{1};
    Y = pos_G_W{2};
    Z = pos_G_W{3};
    U = f_kite_G{1};
    V = f_kite_G{2};
    W = f_kite_G{3};
    quiver3(X,Y,Z,U,V,W) 
end

F_t_W = f_kite_G;

% Loop through all particles
v_w_vec = cell(3,1);
v_w_vec{2} = 0;
v_w_vec{3} = 0;
for idx = 1 : n_t_p
    v_w_vec{1} = v_w_u{idx};
    % Calculate Spring-Damper forces
    tmp = obj.calcSpringDamperForce( c, p_norm{idx}, tether_diff, d, p_d( (idx-1)*3+1:(idx-1)*3+3 ), v_d( (idx-1)*3+1:(idx-1)*3+3) );
    fsd_tot_mat{1,idx} = tmp{1};
    fsd_tot_mat{2,idx} = tmp{2};
    fsd_tot_mat{3,idx} = tmp{3};
    if idx < n_t_p
        fsd_tot_mat(:,idx+1) = obj.calcSpringDamperForce( c, p_norm{idx+1}, tether_diff, d, p_d( idx*3+1:idx*3+3 ), v_d( idx*3+1:idx*3+3) );
    end
    
    % Calculate apparent wind speed of tether particles
    if idx == 1 % Segment that is connected to the ground
        v_a_p = sub_cellMatrix(v_w_vec,mult_cellMatrix({0.5},v_vec(1:3)));
    elseif idx == n_t_p
        v_a_p = sub_cellMatrix(v_w_vec,mult_cellMatrix({0.5},add_cellMatrix(v_vec((idx-2)*3+1:(idx-2)*3+3),v_vec( (idx-1)*3+1:(idx-1)*3+3)))); % take mean value of adjacent particle vel. for segment velocity
    else
        v_a_p = sub_cellMatrix(v_w_vec,mult_cellMatrix({0.5},add_cellMatrix(v_vec((idx-2)*3+1:(idx-2)*3+3),v_vec( (idx-1)*3+1:(idx-1)*3+3)))); % take mean value of adjacent particle vel. for segment velocity
    end
    
    % Calculate segment drag
    D_s_mat(:,idx) = obj.calcSegmentDrag(p_d((idx-1)*3+1:(idx-1)*3+3) , v_a_p, obj.ENVMT.rhos, CD_tether, d_tether); % segment drag
    
    % Calculation of the segment tether force
    % kite_tether_force(:,idx) = F_i
    % D_s_mat = Fa
    % 
    % F_i                      =     Fs,i+1  −   Fs,i−   Fg    +   Fa
    
    if idx < n_t_p
        % Equations of motion for the groud tether particles Newton kite_tether_force = (-f + u -mg )
        if idx+1 == n_t_p
            kite_tether_force(:,idx) = add_cellMatrix(sub_cellMatrix(sub_cellMatrix(fsd_tot_mat(:,idx+1), fsd_tot_mat(:,idx)), Fg),D_s_mat(:,idx));  % Particle dynamics
        else
            if idx == 1
                kite_tether_force(:,idx) = add_cellMatrix(sub_cellMatrix(sub_cellMatrix(fsd_tot_mat(:,2), fsd_tot_mat(:,1)),Fg),D_s_mat(:,1));  % Particle dynamics
            else
                kite_tether_force(:,idx) = add_cellMatrix(sub_cellMatrix(sub_cellMatrix(fsd_tot_mat(:,idx+1) ,fsd_tot_mat(:,idx)),Fg),D_s_mat(:,idx));  % Particle dynamics
            end
        end
    else
        % Flexible tether: last particle of the tether is connected to the
        % aircraft CG.
        max_F_tether =             add_cellMatrix(sub_cellMatrix(sub_cellMatrix(f_kite_G, fsd_tot_mat(:,idx)),Fg),D_s_mat(:,idx)); 
        kite_tether_force(:,idx) = add_cellMatrix(sub_cellMatrix(sub_cellMatrix(f_kite_G, fsd_tot_mat(:,idx)),Fg),D_s_mat(:,idx));  % Particle dynamics
    end
end

ground_tether_force =  fsd_tot_mat(:,1);

end

