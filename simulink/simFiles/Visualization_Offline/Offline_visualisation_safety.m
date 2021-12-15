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
% :Adapted from: Dylan Eijkelhof (d.eijkelhof@tudelft.nl)

function [fig_saf] = Offline_visualisation_safety(safety_last_cycle,Path_last_cycle)
%Plot flight path and safety switching
%
% :param safety_last_cycle: Converged cycle timeseries of safety switching
% :param Path_last_cycle: Converged cycle timeseries of kite position

%------------- BEGIN CODE --------------
%% Flight path & power
fig_saf = figure('Name', 'Flight path, switching coloured');
set(gcf, 'Position',  [100, 100, 1400, 1000])
axes1 = axes('Parent',fig_saf,'Position',[0.13,0.11,0.678584552871244,0.815]);
hold(axes1,'on');
fig_saf.Renderer = 'painters';
view([16 19])
fontsize = 35;

% Flight path coloured by safety switch
sizeScatter = 14;
Safety = safety_last_cycle.Data;
scatter3(-Path_last_cycle.Data(1,:),Path_last_cycle.Data(2,:),-Path_last_cycle.Data(3,:),sizeScatter*ones(size(Safety)),Safety,'filled'); hold on
enhance_plot('Times',fontsize,2,40,0)

%Colorbar
LimitsColorBar = [floor(min(Safety)/1)*1 ceil(max(Safety)/1)*1];

CM = bone((LimitsColorBar(2)-LimitsColorBar(1))*100);
negbar = hot(abs(LimitsColorBar(1))*100+900);
CM(1:abs(LimitsColorBar(1))*100,:) = negbar(1:end-900,:);
CM(abs(LimitsColorBar(1))*100+1:end,:) = flipud(winter(size(CM(abs(LimitsColorBar(1))*100+1:end,:),1)));
colormap(CM)  % Set the colormap of the figure

cb = colorbar('peer',axes1,'Location','eastoutside','FontSize',16);
cb.Label.String = 'Safety [On/Off]';
cb.Label.Interpreter = 'latex';
cb.TickLabelInterpreter = cb.Label.Interpreter;
caxis(LimitsColorBar)
cb.FontSize = fontsize;
cb.Location = 'eastoutside';
cb.Position = [0.874536933823625,0.111493416432584,0.010666666666667,0.815];   
cb.Label.FontSize = fontsize;
cb.Ticks = [0 1];
cb.TickLabels = {'False','True'};


% Graph limits
limitz = [0 500];
limity = [-300 300];
limitx = [-100 600];

xlabel('$$x_W$$ [m]','FontSize',fontsize, 'interpreter', 'latex');
ylabel('$$y_W$$ [m]','FontSize',fontsize, 'interpreter', 'latex');
zlabel('$$z_W$$ [m]','FontSize',fontsize, 'interpreter', 'latex');
set(gca,'TickLabelInterpreter','latex');

X = [limitx(1);limitx(2);limitx(2);limitx(1)];
Y = [limity(1);limity(1);limity(2);limity(2)];
Z = [limitz(1);limitz(1);limitz(1);limitz(1)];
c = [0 1 0; 
0 1 0; 
0 1 0;
0 1 0];
fill3(X,Y,Z, c(:,2),'EdgeColor','none','FaceColor',c(1,:),'FaceAlpha',0.1);

%Path projections on sides/walls of the graph
h(1) = scatter3(-Path_last_cycle.Data(1,:),Path_last_cycle.Data(2,:),zeros(size(-Path_last_cycle.Data(3,:))));
% enhance_plot('Times',fontsize,2,20,0)
h(2) = scatter3(-Path_last_cycle.Data(1,:),zeros(size(Path_last_cycle.Data(2,:)))+limity(2),-Path_last_cycle.Data(3,:));
% enhance_plot('Times',fontsize,2,20,0)
% h(3) = scatter3(zeros(size(Path_last_cycle.Data(2,:)))+limitx(1),Path_last_cycle.Data(2,:),-Path_last_cycle.Data(3,:));
set(h,'MarkerFaceColor',[0.4 0.4 0.4],'MarkerEdgeColor','none','SizeData',.5,'LineWidth',0.01);

grid on
axis equal
zlim(limitz)
ylim(limity)
xlim(limitx)

%------------- END CODE --------------
end

