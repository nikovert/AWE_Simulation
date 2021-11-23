function visualize_vehicle(uu)
    % process input to function
    pn = -uu(1); % inertial North position
    pe = uu(2); % inertial East position
    pd = -uu(3); % inertial down position

    phi = uu(4);
    theta = uu(5);
    psi = uu(6);

    t = uu(7);

    if length(uu) > 7
        Ft = uu(8);
    else
        Ft = 0;
    end

    p_VT_W = [];

    complex_tether_flag = 0;
    persistent tether_handle
    persistent geo_tether_handle
    persistent counter
    persistent p_VT_W_handle

    persistent Vert
    persistent p

    if t==0
        counter = 1;
        scale = 10; 
        [~, V, ~] = rndread('kitemill2.stl');
        Vert = scale*V';
        
        [tether_handle, geo_tether_handle] = drawTether(pn, pe, pd, [], []); hold on
          
        pos_W = [-pn; pe; -pd];
        p_VT_W_handle = drawVirtualTarget(p_VT_W,pos_W,[]);
        
        [p] = drawVehicleBody2(Vert,pn, pe, pd, phi, theta, psi,[],'normal');

        color = [0, 85/255, 85/255] ;
        pathpoints = animatedline('Linestyle','-','color', color, 'Linewidth', 1.5);
        pathpoints2 = animatedline('Marker','.','color', [0 0.5 0], 'Linewidth', 0.2);
        axis([-50 600 -300 300 0 500]);  hold on

        %======== lissajous figure ========
        if 0
            l_tether = norm(uu(1:3));
            [ LemPs ] = updateLissajous( l_tether, LemPsRef );
            Alambda = LemPs.Alambda ;
            Aphi = LemPs.Aphi;
            blambda = LemPs.blambda;
            bphi = LemPs.bphi;
            phi0 = LemPs.phi0;
            theta_vec = 0 : 0.001 : 2*pi;theta_vec = theta_vec';
            L = [Alambda * sin(blambda*theta_vec');
                Aphi    * sin(bphi*theta_vec') + phi0];
            L_W = [cos(L(1,:)).*cos(L(2,:));
                sin(L(1,:)).*cos(L(2,:));
                sin(L(2,:))]*l_tether;
        end
        %=================================
        
        xlabel('$$(m)$$', 'interpreter', 'latex')
        ylabel('$$(m)$$', 'interpreter', 'latex')
        zlabel('$$(m)$$', 'interpreter', 'latex')
        set(gca, 'Fontsize' ,14);
        grid on
        box on;
        axis equal; hold on;
        view(50,30); % set the view angle for figure

        set(gca,'TickLabelInterpreter','latex');
        
        c = colorbar('southoutside');
        c.Label.Interpreter = 'latex';
        c.TickLabelInterpreter = c.Label.Interpreter;
        c.Limits = [0 1];
        c.Label.String = 'Tether Force ($\%$ of  $F_{t,\max}$)';
        set(c, 'Ticks',[0,0.25,0.5,0.75,1]);
        set(c, 'TickLabels',{'$0\%$','$25\%$','$50\%$','$75\%$','Tether Rupture'});
    else
        axis([-50 600 -300 300 0 500]);  hold on
        drawVehicleBody2(Vert,pn, pe, pd, phi, theta, psi,p);
        
        Ft_max = 1861;
        if Ft < Ft_max
            force_color = Ft/Ft_max;
            scatter3(-pn, pe, -pd,14, force_color, 'filled');
        else
            scatter3(-pn, pe, -pd,40, 'red', '*');
        end
        cstm_plot3(-pn, pe, -pd, 'MarkerSize', 5, 'Color', 'red', 'LineWidth', 3);
        
        [tether_handle, geo_tether_handle] = drawTether(pn, pe, pd, tether_handle, geo_tether_handle);

        pos_W = [-pn;pe;-pd];
        if norm(pos_W) > 1e-3
            p_VT_W_handle = [];drawVirtualTarget(p_VT_W,pos_W, p_VT_W_handle);
        else
            p_VT_W_handle = [];
        end
        if 1
            deltaC = 1;
            if mod(counter-1, deltaC) == 0
                if counter == 1
                    print([eval('pwd'),'/video/','vidPic_',num2str(counter)], '-dpng', '-r300');
                else
                    print([eval('pwd'),'/video/','vidPic_',num2str((counter-1)/deltaC)], '-dpng', '-r300');
                end
            end
        end

        %======== lissajous figure ========
        if 0
            l_tether = norm(uu(1:3));
            [ LemPs ] = updateLissajous( l_tether, LemPsRef );
            Alambda = LemPs.Alambda ;
            Aphi = LemPs.Aphi;
            blambda = LemPs.blambda;
            bphi = LemPs.bphi;
            phi0 = LemPs.phi0;
            theta_vec = 0 : 0.001 : 2*pi;theta_vec = theta_vec';
            L = [Alambda * sin(blambda*theta_vec');
                Aphi    * sin(bphi*theta_vec') + phi0];
            L_W = [cos(L(1,:)).*cos(L(2,:));
                sin(L(1,:)).*cos(L(2,:));
                sin(L(2,:))]*l_tether;
        end
        %=================================
        counter = counter + 1;

    end
end

function handle = cstm_plot3(X,Y,Z, varargin)
    handle = [];
    return;

    persistent wgs84
    persistent g
    gndlat = 51.779438;
    gndlon = -1.287673;
    gndalt = 60;
    if isempty(g)
        fig = uifigure;
        g = geoglobe(fig);
        hold(g,"on")
        
        % Plot groundstation
        geoplot3(g,gndlat,gndlon,gndalt,"co", ...
            "LineWidth",6, ...
            "MarkerSize",1)
        
        campos(g,51.7763,-1.2964,500)
        camheading(g,70)
        campitch(g,-12)
    end
    
    if isempty(wgs84)
        wgs84 = wgs84Ellipsoid;
    end

    [tlat,tlon,tht] = enu2geodetic(X,Y,Z,gndlat,gndlon,gndalt,wgs84);
    
    if length(tlat) == 2
        linarg = '-';
    else
        linarg = 'o';
    end
    
    if nargin < 4
        handle = geoplot3(g,tlat,tlon,tht,linarg, "HeightReference","ellipsoid");
    else
        handle = geoplot3(g,tlat,tlon,tht, linarg, ...
            "HeightReference","ellipsoid", varargin{:});
    end
end

function set_geo(handle, X, Y, Z)
    persistent wgs84
    gndlat = 51.779438;
    gndlon = -1.287673;
    gndalt = 60;
    if isempty(wgs84)
        wgs84 = wgs84Ellipsoid;
    end

    [tlat,tlon,tht] = enu2geodetic(X,Y,Z,gndlat,gndlon,gndalt,wgs84);
    set(handle,'LatitudeData',tlat,'LongitudeData',tlon,'HeightData',tht);
end
function handleParticles = drawParticleTether(p_t_x,p_t_y,p_t_z, handleParticles)
    if isempty(handleParticles)
        handleParticles = cstm_plot3( p_t_x, p_t_y,p_t_z , '-o', 'color', [0.1 0.1 0.1], 'Markersize',2, 'Linewidth', 1, 'MarkerFaceColor', [0.3 0.3 0.3]); hold on;
    else
        set_geo(handleParticles, p_t_x,p_t_y,p_t_z)
        %set(handleParticles,'XData',p_t_x,'YData',p_t_y,'ZData',p_t_z);
    end
    drawnow;
end

function handleTempPath = drawTempPath(x,y,z, handleTempPath)
    if isempty(handleTempPath)
        color =  [   128, 128, 128   ]/255;
        handleTempPath = cstm_plot3( x, y,z , '-', 'color', color,'Linewidth', 2); hold on;
    else
        set_geo(handleTempPath, x,y,z)
        %set(handleTempPath,'XData',x,'YData',y,'ZData',z);
    end
    drawnow;
end

function referencePath = drawReferencePath(x,y,z, referencePath)
    if isempty(referencePath)
        referencePath = cstm_plot3( x, y,z , '-', 'color', [0.3 0.3 0.3],'Linewidth', 1.5); hold on;
    else
        set_geo(referencePath, x,y,z)
        %set(referencePath,'XData',x,'YData',y,'ZData',z);
    end
    drawnow;
end

function [handleT, handleT_geo]  = drawTether(pn, pe, pd, handleT, handleT_geo)
    M_EO = [-1, 0, 0; 0, 1, 0; 0, 0, -1];
    pts = M_EO * [pn;pe;pd];
    
    if isempty(handleT)
        handleT = plot3([0 pts(1)], [0 pts(2)], [0 pts(3)], 'Color', 'black');
    else
        set(handleT,'XData',[0 pts(1)],'YData',[0 pts(2)],'ZData',[0 pts(3)] );
    end
    
    if isempty(handleT_geo)
        handleT_geo = cstm_plot3([0 pts(1)], [0 pts(2)], [0 pts(3)], 'Color', 'black');
    else
        set_geo(handleT_geo, [0 pts(1)], [0 pts(2)], [0 pts(3)]);
    end
    drawnow;
end

function p_VT_W_handle = drawVirtualTarget(p_VT_W,pos_W, p_VT_W_handle)
    if isempty(p_VT_W_handle)
        p_VT_W_handle = [];%plot3(p_VT_W(1), p_VT_W(2), p_VT_W(3), '+b');
        %quiver3( pos_W(1), pos_W(2), pos_W(3), p_VT_W(1)*50, p_VT_W(2)*50,p_VT_W(3)*50,'b')%
    else
        set(p_VT_W_handle,'XData',p_VT_W(1),'YData',p_VT_W(2),'ZData',p_VT_W(3));
    end
    drawnow;
end

function [p] = drawVehicleBody2(V,pn, pe, pd, phi, theta, psi,p, mode)
    %V = rotate( V, phi+pi, theta, psi+pi );
    V = rotate( V, -phi-pi, -theta, psi+pi );
    % Stabilizer
    %V_stab = rotate(V_stab, phi, theta, psi); % body frame into NED frame
    V = translate(V, pn , pe, pd);
    M_EO = [-1, 0, 0; 0, 1, 0; 0, 0, -1];
    V = M_EO*V; % Transform from NED into E frame

    if isempty(p)
        [F, ~, C] = rndread('kitemill2.stl');
        col =  [ 57/255, 106/255, 177/255  ];
        p = patch('faces', F, 'vertices' ,V');
        set(p, 'facec', 'b');              % Set the face color (force it)
       % set(p, 'facec', 'flat');            % Set the face color flat
        set(p, 'FaceVertexCData', C);       % Set the color (from file)
        %set(p, 'facealpha',.4)             % Use for transparency
        set(p, 'EdgeColor','none');
        set(p, 'FaceLighting', 'none'); 
        set(p, 'FaceColor', col ); 
       % light  
        %('Position',[0 0 1]);% add a default light
        view(45,45);
        daspect([1 1 1])   
    else
        %set(p, 'Xdata', xstab, 'Ydata', ystab, 'Zdata', zstab);
        set(p, 'Vertices', V');
        %set(h_wing, 'Xdata', xwing, 'Ydata', ywing, 'Zdata', zwing);
        %set(h_plane, 'Xdata', xplane, 'Ydata', yplane, 'Zdata', zplane);
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

function [V_stab,V_wing,V_plane] = defineVehicleBody
    scale =5;
    xplane = [];
    yplane = [];
    zplane = [];
    R = [ acos([0;0.4;0.8;1]) ];
    [x,y,z] = cylinder(R,10);
    z(end,:) = [3.*z(end,:) + 4.*z(end-1,:)]./7;
    aux=x;x=-z;z=aux;
    x = x./max(max(abs(x)));
    y = y./max(max(abs(y)));
    z = z./max(max(abs(z)));
    x = flipud(x);
    y = flipud(y);
    z = flipud(z);
    xnose = x.*2.50 + 3.00;
    ynose = y.*0.80 + 0.00;
    znose = z.*0.85 + 3.15;
    xplane = [ xplane; xnose ];
    yplane = [ yplane; ynose ];
    zplane = [ zplane; znose ];
    xfus1 = x(end,:).*0.00 + 8.5;
    yfus1 = ynose(end,:);
    zfus1 = znose(end,:);
    xplane = [ xplane; xfus1 ];
    yplane = [ yplane; yfus1 ];
    zplane = [ zplane; zfus1 ];
    xfus2 = x(end,:).*0.00 + 12.5;
    yfus2 = y(end,:).*0.30;
    zfus2 = z(end,:).*0.55 + 3.7;
    xplane = [ xplane; xfus2 ];
    yplane = [ yplane; yfus2 ];
    zplane = [ zplane; zfus2 ];
    xtail = [13.8; 14.1; 14.8]*ones(size(x(end,:)));
    ytail = [y(end,:).*0.10 + 0.00; y(end,:).*0.10 + 0.00; y(end,:).*0.00 + 0.00];
    ztail = [z(end,:).*1.45 + 4.85; z(end,:).*1.40 + 4.90; z(end,:).*0.00 + 6.30];
    xplane = [ xplane; xtail ];
    yplane = [ yplane; ytail ];
    zplane = [ zplane; ztail ];
    % xwing = [ ...
    %     7.0   7.25 8.0  7.25 7.0; ...
    %     6.0   7.0  9.5 7.0  6.0; ...
    %     7.0   7.25 8.0  7.25 7.0; ...
    %     ];
    xwing = [ ...
        7.0   7.25 8.0  7.25 7.0; ...
        6.5   7.0  9 7.0  6.5; ...
        7.0   7.25 8.0  7.25 7.0; ...
        ];
    ywing = 1.5*[ ...
        -8.0 -8.0 -8.0 -8.0 -8.0; ...
        0.0   0.0  0.0  0.0  0.0; ...
        8.0   8.0  8.0  8.0  8.0];
    zwing = 0.2 + [ ...
        3.0   3.0  3.0  3.0  3.0; ...
        2.6   2.55 2.4  2.45 2.6; ...
        3.0   3.0  3.0  3.0  3.0; ...
        ];
    xstab = (xwing - 7.0).*(1.5/3.5) + 13;
    ystab = (ywing - 0.0).*(3.0/8.0) +  0;
    zstab = [ ...
        4.5   4.5  4.5  4.5  4.5; ...
        4.1   4.1  4.1  4.1  4.1; ...
        4.5   4.5  4.5  4.5  4.5; ...
        ];

    % cg correction
    x_cg   = 7;
    xstab  = xstab  - x_cg;
    xwing  = xwing  - x_cg;
    xplane = xplane - x_cg;
    y_cg   = 0;
    ystab  = ystab  - y_cg;
    ywing  = ywing  - y_cg;
    yplane = yplane - y_cg;
    z_cg   = 3.2;
    zstab  = zstab  - z_cg;
    zwing  = zwing  - z_cg;
    zplane = zplane - z_cg;

    V_stab = scale*[ reshape(-xstab, 1, 15);reshape(ystab, 1, 15);reshape(-zstab, 1, 15)];
    V_wing = scale*[ reshape(-xwing, 1, 15);reshape(ywing, 1, 15);reshape(-zwing, 1, 15)];
    V_plane = scale*[ reshape(-xplane, 1, 99);reshape(yplane, 1, 99);reshape(-zplane, 1, 99)];
end
