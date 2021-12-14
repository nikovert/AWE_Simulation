% Copyright 2021 Delft University of Technology
%
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
%
%      http://www.apache.org/licenses/LICENSE-2.0
%
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.

function [figF, fig_force] = Offline_visualisation_force(Force_last_cycle,Path_last_cycle, Ft_max, rupture)
%Plot power curve and flight path
%
% :param Force_last_cycle: Converged cycle timeseries of power
% :param Path_last_cycle: Converged cycle timeseries of kite position
% :param Ft_max: max alowed tether force
% :param logic if rupture occured
% :returns: None
%
% | Other m-files required: enhance_plot.m
% | Subfunctions: none
% | MAT-files required: none
%
% :Revision: 01-September-2020
% :Author: Dylan Eijkelhof (d.eijkelhof@tudelft.nl)

%------------- BEGIN CODE --------------
if nargin < 4
    rupture = false;
end
%% Power
figF = figure('Name', 'Tether Force');
Force = Force_last_cycle.Data(:,2)./1e3;
PTime = Force_last_cycle.Time-Force_last_cycle.Time(1);

% Level for Color Change
lev = Ft_max./1e3;
% Find points above the level
aboveLine = (Force>lev);
% Create 2 copies of y
validFore = Force;

% Set the values you don't want to get drawn to nan
validFore(aboveLine) = NaN;


plot(PTime,validFore,'Color','#77AC30','Marker','none','MarkerIndices',1:100:length(validFore))
hold on; box on; grid on; axis tight;
%plot(PTime,rupture,'Color','#A2142F','Marker','none','MarkerIndices',1:100:length(rupture))

plot([0 PTime(end)],[mean(Force) mean(Force)],'--','Color','#EDB120')

enhance_plot('times',16,1.5,4,3)
% legend('Power','', 'Average power','Peak power','Location','best')
if mean(Force) < 1
    mf = [num2str(mean(Force)*1e3,3) 'N'];
else
    mf = [num2str(mean(Force),3) ' kN'];
end
legend('Tether Force', ['Average Force (' mf ')'],'Location','best')
ylabel('Tether Force [kN]')
xlabel('Time [s]')
%     saveas(figPO,'power_cycle','epsc')

%% Flight path & power
fig_force = figure('Name', 'Flight path, tether force' );
set(gcf, 'Position',  [100, 100, 1400, 1000])
axes1 = axes('Parent',fig_force,'Position',[0.13,0.11,0.678584552871244,0.815]);
hold(axes1,'on');
fig_force.Renderer = 'painters';
view([16 19])
% view([138 24])
% set(gca,'XColor', 'none','YColor','none','ZColor','none')
% set(gcf, 'color', 'white');
% set(gca, 'color', 'white');
fontsize = 35;
% Flight path coloured by power
sizeScatter = 14;
Force = Force_last_cycle.Data(:,2)./Ft_max;
scatter3(-Path_last_cycle.Data(1,:),Path_last_cycle.Data(2,:),-Path_last_cycle.Data(3,:),sizeScatter*ones(size(Force)),Force,'filled'); hold on
enhance_plot('Times',fontsize,2,40,0)

if rupture
    % Draw tether rupture point
    scatter3(-Path_last_cycle.Data(1,end),Path_last_cycle.Data(2,end),-Path_last_cycle.Data(3,end), 200, 'r*');
end

% Colorbar
LimitsColorBar = [0 1];
 
CM = bone((LimitsColorBar(2)-LimitsColorBar(1))*100);
negbar = hot(abs(LimitsColorBar(1))*100+900);
CM(1:abs(LimitsColorBar(1))*100,:) = negbar(1:end-900,:);
CM(abs(LimitsColorBar(1))*100+1:end,:) = turbo(size(CM(abs(LimitsColorBar(1))*100+1:end,:),1));
colormap(CM)  % Set the colormap of the figure

cb = colorbar('peer',axes1,'Location','eastoutside','FontSize',16);
cb.Label.Interpreter = 'latex';
cb.Label.String = 'Force [\% of $$F_{\mathrm{rupture}}$$]';
cb.TickLabelInterpreter = cb.Label.Interpreter;
caxis(LimitsColorBar)
cb.FontSize = fontsize;
cb.Location = 'eastoutside';
cb.Position = [0.874536933823625,0.111493416432584,0.010666666666667,0.815];   
cb.Label.FontSize = fontsize;

cb.Ticks = [0 0.5 1];
cb.TickLabels = {'0\%','50\%', '100\%'};

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

% X = [limitx(1);limitx(1);limitx(1);limitx(1)];
% Y = [limity(1);limity(1);limity(2);limity(2)];
% Z = [limitz(1);limitz(2);limitz(2);limitz(1)];
% c = [200/255 200/255 200/255;
% 0 0 0;
% 0 0 0;
% 0 0 0];
% fill3(X,Y,Z, c(:,2),'EdgeColor','none','FaceColor',c(1,:),'FaceAlpha',0.1);
% 
% X = [limitx(1);limitx(2);limitx(2);limitx(1)];
% Y = [limity(2);limity(2);limity(2);limity(2)];
% Z = [limitz(1);limitz(1);limitz(2);limitz(2)];
% fill3(X,Y,Z, c(:,2),'EdgeColor','none','FaceColor',c(1,:),'FaceAlpha',0.1);

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

