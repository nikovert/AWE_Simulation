function visualize(x, u, windDirection)
%VISUALIZE Summary of this function goes here
%   Detailed explanation goes here
    persistent Vertices
    persistent Faces
    persistent facecolors
    persistent vehicle_handle
    persistent pathpoints
    
    if nargin < 3
        windDirection = pi;
    end
    
    long    = x(1);
    lat     = x(2);
    h_tau   = x(3);
    va      = x(4);
    chi_a   = x(5);
    gamma_a = x(6);
    
    alpha   = u(1);
    beta    = 0;
    mu_a    = u(2);
      
    M_WO = [cos(windDirection), sin(windDirection), 0;
        sin(windDirection), -cos(windDirection), 0;
        0, 0, -1];
    
    M_AbarO = [cos(chi_a) .* cos(gamma_a), sin(chi_a) .* cos(gamma_a), -sin(gamma_a); 
              -sin(chi_a),                 cos(chi_a),                             0; 
               cos(chi_a) .* sin(gamma_a), sin(chi_a) .* sin(gamma_a), cos(gamma_a)];
           
    M_OAbar = M_AbarO';
    M_AbarA = [1,0,0;
               0, cos(mu_a), -sin(mu_a);
               0, sin(mu_a), cos(mu_a)];

    M_AB = [cos(alpha) .* cos(beta), sin(beta), sin(alpha) .* cos(beta); 
        -cos(alpha) .* sin(beta), cos(beta), -sin(alpha) .* sin(beta); 
        -sin(alpha), 0, cos(alpha)];

    M_OB = M_OAbar * M_AbarA * M_AB;
    
    phi =  atan2( M_OB(3,2), M_OB(3,3) ); 
    theta = -asin( M_OB(3,1) ); 
    psi = atan2( M_OB(2,1), M_OB(1,1) ); 
    
    [pos_x,pos_y,pos_z] = sph2cart(long,lat,h_tau);
    pos_O = [pos_x;pos_y;pos_z];
    
    [Vertices, Faces, facecolors] = defineVehicleBody;
    
    vehicle_handle = drawVehicleBody(Vertices, Faces, facecolors,...
         -pos_O(1), pos_O(2), -pos_O(3), phi, theta, psi, ...
        vehicle_handle, windDirection);
    pathpoints = animatedline('color', [0.5 0.5 0.5], 'Linewidth', 1);
    addpoints(pathpoints, pos_O(1), pos_O(2), pos_O(3));
    drawnow update 
end

function handle = drawVehicleBody(V,F,patchcolors,pn, pe, pd, phi, theta, psi,handle, windDirection)
    if nargin < 11
        windDirection = pi;
    end
    V = rotate(V, phi, theta, psi); % body frame into NED frame
    V = translate(V, pn , pe, pd);
    
    M_WO = [cos(windDirection), sin(windDirection), 0;
        sin(windDirection), -cos(windDirection), 0;
        0, 0, -1];
    V = M_WO*V;
    if isempty(handle) || ~isvalid(handle)
        handle = patch('Vertices', V', 'Faces', F, ....
            'FaceVertexCData', patchcolors,...
            'FaceColor', 'flat');%,...
        %'EraseMode', mode);
        %axis([-50 300 -300 300 -100 300]); hold on;
        view(30,60)
    else
        set(handle,'Vertices', V', 'Faces', F);
        drawnow
    end
end

function pts = rotate(pts, phi, theta, psi)

    % From B 2 O
    pts = [ cos(psi)*cos(theta), cos(psi)*sin(phi)*sin(theta) - cos(phi)*sin(psi), sin(phi)*sin(psi) + cos(phi)*cos(psi)*sin(theta);
        cos(theta)*sin(psi), cos(phi)*cos(psi) + sin(phi)*sin(psi)*sin(theta), cos(phi)*sin(psi)*sin(theta) - cos(psi)*sin(phi);
        -sin(theta),                              cos(theta)*sin(phi),                              cos(phi)*cos(theta)] * pts;
end

function pts = translate(pts, pn, pe, pd)
    pts = pts + repmat([pn;pe;pd], 1, size(pts,2));
end

function [V,F,facecolors] = defineVehicleBody
    scale = 15;
    b_wing = 3.7*scale;
    c_wing = 0.22*scale;
    L_fuse = 1*scale;
    b_fuse = 0.1*scale;
    c_emp = 0.1*scale;
    b_emp = 0.4*scale;
    h_rud = 0.3*scale;
    h_fuse = 0.1*scale;
    x_nose = 0.5*scale;
    x_tail = 1.5*scale;
    b_fuse = 0.1*scale;
    xp = 0.1*scale;
    x_LE = 0.1*scale;
    c_rud = 0.1*scale;
    V = [...
        x_nose, 0, 0; %1
        x_nose-xp, b_fuse/2, -h_fuse/2; %2
        x_nose-xp, b_fuse/2, +h_fuse/2; %3
        x_nose-xp, -b_fuse/2, +h_fuse/2; %4
        x_nose-xp, -b_fuse/2, -h_fuse/2; %5
        -x_tail, 0, 0; %6
        x_LE, b_wing/2, 0;%7
        x_LE-c_wing, b_wing/2, 0;%8
        x_LE-c_wing, -b_wing/2, 0; %9
        x_LE, -b_wing/2, 0;  %10
        -x_tail, -b_emp/2, 0;... % 11
        -(x_tail+c_emp), -b_emp/2, 0;...% 12
        -(x_tail+c_emp), b_emp/2, 0;...% 13
        -x_tail, b_emp/2, 0;...% 14
        -(x_tail+c_emp), 0, 0;... % 15
        -(x_tail+c_emp+c_rud), 0, 0;...% 16
        -(x_tail+c_emp+c_rud), 0, -h_rud;...% 17
        -(x_tail+c_emp), 0, -h_rud;...% 18
        ]';
    F = [...
        1,2,3,3;... % pyramid
        3,4,1,1;...
        1,4,5,5;...
        1,2,5,5;...
        3,4,6,6;...
        4,5,6,6;...
        2,5,6,6;...
        2,3,6,6;...
        2,5,6,6;...
        7,8,9,10;...
        11,12,13,14;...
        15,16,17,18;...
        ];
    facecolors = [1,0,0];
end