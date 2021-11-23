function [x_dotdot_TetherG2K,l_s_dot, ground_tether_force]    = pure_tether_model_function(x_dot, x,  v_w,...
    d_tether, CD_tether, v_ro_set,l_s, c0, d0, rho_t, g_vec, rho_air,...
    x_kite_A,x_kite_C,x_kite_D ,xdot_kite_A,xdot_kite_C,xdot_kite_D, ls_AKcu, ls_CKcu,  ls_DKcu, m_kcu, landedFlag   )
%% Readme
% Equations are adapted from Fechner et al, Dynamic Model of a Pumping Kite
% Power System, Renewable Energy, 2015.
% Implementation: Sebastian Rapp, Wind Energy Institute, Faculty of
% Aerospace Engineering, TU Delft
% Mail: s.rapp@tudelft.nl
% Last change: 19.01.2018
%
% General description:
% This functions calculates the accelerations of the tether particle
% segment based on the acting forces at each particle.
%
% Inputs:
% x_dot: Tether particle velocities (world frame i.e. inertial frame)
% x: Tether particle position (world frame)
% v_w: wind vector (world frame)
% d_tether: tether diameter
% CD_tether: tether drag coefficient
% v_ro_set: reelout speed of the tether
% l_s: tether segment length
% c0: intial spring constant of a tether segment
% d0: initial damping constant of a tether segment
% rho_t: tether mass density (i.e. mass per m)
% g_vec: gravity vector
% rho_air: air density
% x_kite_A,C,D: position bridle attachment points
% xdot_kite_A,C,D: velocities of bridle attachment points
% ls_AKcu,CKcu,DKcu: bridle lengths
% m_kcu: kcu weight
%
% Outputs:
% x_dotdot_TetherG2K: accelertions of the tether particles
% l_s_dot: tether length increases with this speed
% ground_tether_force: tether force as measured at the ground

%#codegen

n_t_p = 5; % amount of particles
l_s_dot = v_ro_set/n_t_p; % Particle speeds

% Adapt stiffness and damping with respect to current segment length
c = c0/l_s;
d = d0/l_s;

% Mass of tether segments with  respect to current segment length
m = l_s * rho_t;
m_vec = m * ones(n_t_p,1);
m_vec(n_t_p) = m_kcu + (l_s/2 * rho_t);

p_vec = x;
v_vec = x_dot;
p_norm = zeros(n_t_p,1);

% Vector with distances between particles
p_d = p_vec - [0;0;0; p_vec(1:3*n_t_p-3)];

% Vector with velocity differences
v_d = v_vec - [0;0;0; v_vec(1:3*n_t_p-3)];

% calculate absolute distances between points
for iter = 1 : n_t_p
    p_norm(iter) = norm( p_d( (iter-1)* 3 + 1 : (iter-1) * 3 + 3 ) );
end

%%------------------------------- Set up equations of motion-------------------------------
% Initialize matrices
fsd_tot_mat = zeros(3,n_t_p);
x_dotdot = zeros(3,n_t_p);
D_s_mat = zeros(3,n_t_p);

%% Set up equations of motions of ground to kite tether

% Boundary tether force from kite to particle below kite (in the 4p kite model the kcu is connected to the
% points A, C and D

f_kite_A = calcSpringDamperForce( c0,norm( x_kite_A - p_vec(3*(n_t_p) -2 : 3*(n_t_p) ) ), ls_AKcu, d0,...
    x_kite_A - p_vec(3*(n_t_p) -2 : 3*(n_t_p) ) ,...
    xdot_kite_A - v_vec( 3*(n_t_p) -2 : 3*(n_t_p)  ), 1 );
f_kite_C= calcSpringDamperForce( c0,norm( x_kite_C - p_vec(3*(n_t_p) -2 : 3*(n_t_p) ) ), ls_CKcu, d0,...
    x_kite_C - p_vec(3*(n_t_p) -2 : 3*(n_t_p) ) ,...
    xdot_kite_C - v_vec( 3*(n_t_p) -2 : 3*(n_t_p) ), 1 );
f_kite_D = calcSpringDamperForce( c0,norm( x_kite_D - p_vec(3*(n_t_p) -2 : 3*(n_t_p) ) ), ls_DKcu, d0,...
    x_kite_D - p_vec(3*(n_t_p) -2 : 3*(n_t_p) ) ,...
    xdot_kite_D - v_vec( 3*(n_t_p) -2 : 3*(n_t_p) ), 1 );

v_a_p_save = zeros(3, n_t_p);

% Loop through all particles
for idx = 1 : n_t_p
    v_w_vec = v_w( idx,: )';
    % Calculate Spring-Damper forces
    fsd_tot_mat(:,idx) = calcSpringDamperForce( c, p_norm(idx), l_s, d, p_d( (idx-1)*3+1:(idx-1)*3+3 ), v_d( (idx-1)*3+1:(idx-1)*3+3), idx );
    if idx < n_t_p
        fsd_tot_mat(:,idx+1) = calcSpringDamperForce( c, p_norm(idx+1), l_s, d, p_d( idx*3+1:idx*3+3 ), v_d( idx*3+1:idx*3+3), idx );
    end
    
    % Calculate apparent wind speed of tether particles
    if idx == 1 % Segment that is connected to the ground
        v_a_p = v_w_vec - 0.5*x_dot(1:3) ;
    elseif idx == n_t_p
        v_a_p = v_w_vec - 0.5*( x_dot((idx-2)*3+1:(idx-2)*3+3) + x_dot( (idx-1)*3+1:(idx-1)*3+3) ); % take mean value of adjacent particle vel. for segment velocity
    else
        v_a_p = v_w_vec - 0.5*( x_dot((idx-2)*3+1:(idx-2)*3+3) + x_dot( (idx-1)*3+1:(idx-1)*3+3) ); % take mean value of adjacent particle vel. for segment velocity
    end
    v_a_p_save(1:3, idx) = v_a_p;
    
    % Calculate segment drag
    D_s_mat(:,idx) = calcSegmentDrag(p_d((idx-1)*3+1:(idx-1)*3+3) , v_a_p, rho_air, CD_tether, d_tether); % segment drag
    if idx < n_t_p
        % Equations of motion for the groud tether particles Newton x_dotdot = 1/m * (-f + u -mg )
        if idx+1 == n_t_p
            x_dotdot(:,idx) = 1/m_vec(idx) * ( -fsd_tot_mat(:,idx) + 1*fsd_tot_mat(:,idx+1) - m_vec(idx)*g_vec  + D_s_mat(:,idx) );  % Particle dynamics
        else
            if idx == 1
                x_dotdot(:,idx) = 1/m_vec(idx) * ( -fsd_tot_mat(:,1) + fsd_tot_mat(:,2) - m_vec(1)*g_vec  + D_s_mat(:,1) );  % Particle dynamics
            else
                x_dotdot(:,idx) = 1/m_vec(idx) * ( -fsd_tot_mat(:,idx) + 1*fsd_tot_mat(:,idx+1) - m_vec(idx)*g_vec  + D_s_mat(:,idx) );  % Particle dynamics
            end
        end
    else
        % Flexible tether
        x_dotdot(:,idx) = 1/m_vec(idx) * ( -1*fsd_tot_mat(:,idx) + ...
            1*f_kite_A + 1*f_kite_C + 1*f_kite_D - m_vec(idx)*g_vec  + D_s_mat(:,idx) );  % Particle dynamics
    end
    
    ground_tether_force = norm( fsd_tot_mat(:,1) );
    
end

x_dotdot_TetherG2K = reshape(x_dotdot, [3*n_t_p,1]);
if landedFlag
    x_dotdot_TetherG2K = 0 * x_dotdot_TetherG2K;
end

end
