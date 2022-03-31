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

%% Load workspace
addpath('../')
addToolbox;
load('main_run.mat')

%% Plot Value function
if isfield(HJIextraArgs, 'keepLast') && HJIextraArgs.keepLast
    HJIextraArgs.data0 = data0;
end
plotHJIPDE(grid, data, HJIextraArgs)

%% Functions
function extraOuts = plotHJIPDE(g, data, extraArgs)
    clns = repmat({':'}, 1, g.dim);
    if isfield(extraArgs, 'keepLast') && extraArgs.keepLast
        [gPlot, dataPlot] = proj(g, extraArgs.data0, ~extraArgs.visualize.plotData.plotDims, extraArgs.visualize.plotData.projpt);
    else
        [gPlot, dataPlot] = proj(g, data(clns{:}, 1), ~extraArgs.visualize.plotData.plotDims, extraArgs.visualize.plotData.projpt);
    end
    pDims = nnz(extraArgs.visualize.plotData.plotDims);
    
    % Create figure
    figure1 = figure('Color',[1 1 1]);

    % Create axes
    axes1 = axes('Parent',figure1);
    hold(axes1,'on');
    
    % Set level set slice
    if isfield(extraArgs.visualize, 'sliceLevel')
        sliceLevel = extraArgs.visualize.sliceLevel;
    else
        sliceLevel = 0;
    end
    
    % Set defaults
    eAT_visSetIm.sliceDim = gPlot.dim;
    eAT_visSetIm.applyLight = false;
    if isfield(extraArgs.visualize, 'lineWidth')
        eAT_visSetIm.LineWidth = extraArgs.visualize.lineWidth;
        eAO_visSetIm.LineWidth = extraArgs.visualize.lineWidth;
    else
        eAO_visSetIm.LineWidth = 2;
    end
    
    view3D = 0;
    needLight = true;
    
    %---Visualize Initial Value Set----------------------------------------
    if isfield(extraArgs.visualize, 'initialValueSet') &&...
            extraArgs.visualize.initialValueSet
        
        % If we're making a 3D plot, mark so we know to view this at an
        % angle appropriate for 3D
        if gPlot.dim >= 2
            view3D = 1;
        end
        
        % Set up default parameters
        if ~isfield(extraArgs.visualize,'plotColorVS0')
            extraArgs.visualize.plotColorVS0 = 'g';
        end
        
        extraOuts.hVS0 = visSetIm(...
            gPlot, dataPlot, extraArgs.visualize.plotColorVS0,...
            sliceLevel, eAT_visSetIm);
        
        set(extraOuts.hVS0, 'DisplayName','V(x,0)','LineWidth',1,...
            'LineColor',[0 1 0],...
            'LevelList',0)
        
        if isfield(extraArgs.visualize,'plotAlphaVS0')
            extraOuts.hVS0.FaceAlpha = extraArgs.visualize.plotAlphaVS0;
        end
    end
    
%---Visualize Initial Value Function-----------------------------------
    if isfield(extraArgs.visualize, 'initialValueFunction') &&...
            extraArgs.visualize.initialValueFunction
        
        % Set up default parameters
        if ~isfield(extraArgs.visualize,'plotColorVF0')
            extraArgs.visualize.plotColorVF0 = 'g';
        end
        
        if ~isfield(extraArgs.visualize,'plotAlphaVF0')
            extraArgs.visualize.plotAlphaVF0 = .5;
        end
        
        % Visualize Initial Value function (hVF0)
        [extraOuts.hVF0]= visFuncIm(gPlot,dataPlot,...
            extraArgs.visualize.plotColorVF0,...
            extraArgs.visualize.plotAlphaVF0);
        
        set(extraOuts.hVF0, 'DisplayName','zero level-set of V(x,0)',...
            'Parent',axes1,...
            'FaceLighting','phong',...
            'FaceAlpha',0.5,...
            'FaceColor',[0 1 0],...
            'EdgeColor','none')
    end
    
    [gPlot, dataPlot] = proj(g, data(clns{:}, size(data,g.dim+1)), ~extraArgs.visualize.plotData.plotDims, extraArgs.visualize.plotData.projpt);
    
    
    %---Visualize Value Set------------------------------------------------
    if isfield(extraArgs.visualize, 'valueSet') &&...
            extraArgs.visualize.valueSet
        
        if ~isfield(extraArgs.visualize,'plotColorVS')
            extraArgs.visualize.plotColorVS = 'b';
        end
        
        extraOuts.hVS = visSetIm(gPlot, dataPlot, ...
            extraArgs.visualize.plotColorVS, sliceLevel, eAT_visSetIm);
        
        set(extraOuts.hVS, 'DisplayName','zero level-set of V(x,-0.1)',...
            'LineWidth',1,...
            'LineColor',[0 0 1],...
            'LevelList',0)
    end
    
    %---Visualize Value Function-------------------------------------------
    if isfield(extraArgs.visualize, 'valueFunction') && ...
            extraArgs.visualize.valueFunction
        % If we're making a 3D plot, mark so we know to view this at an
        % angle appropriate for 3D
        if gPlot.dim >= 2
            view3D = 1;
        end
        
        % Set up default parameters
        if ~isfield(extraArgs.visualize,'plotColorVF')
            extraArgs.visualize.plotColorVF = 'b';
        end
        
        if ~isfield(extraArgs.visualize,'plotAlphaVF')
            extraArgs.visualize.plotAlphaVF = .5;
        end
        
        % Visualize Value function (hVF)
        [extraOuts.hVF]= visFuncIm(gPlot,dataPlot,...
            extraArgs.visualize.plotColorVF,...
            extraArgs.visualize.plotAlphaVF);
        
        set(extraOuts.hVF, 'DisplayName','V(x,-0.1)','Parent',axes1,...
            'FaceLighting','phong',...
            'FaceAlpha',0.5,...
            'FaceColor',[0 0 1],...
            'EdgeColor','none')
        
    end
    %% General Visualization Stuff
    
    % Create light
    light('Parent',axes1,'Position',[-30 -30 -30],'Style','local');

    % Create zlabel
    zlabel({'V(x,t)'},'Interpreter','latex');

    % Create ylabel
    ylabel('$\sigma$','Interpreter','latex');

    % Create xlabel
    xlabel({'s'},'Interpreter','latex');

    % Create title
    title('t = -0.1 s','Interpreter','latex');

    view(axes1,[-147.641667373427 14.2061313974064]);
    axis(axes1,'tight');
    axis(axes1,'square');
    hold(axes1,'off');
    % Set the remaining axes properties
    set(axes1,'FontSize',15,'GridAlpha',0.6,'GridLineStyle','--',...
        'MinorGridAlpha',0.2,'TickLabelInterpreter','latex','XMinorGrid','on',...
        'XTick',...
        [0 1.5707963267949 3.14159265358979 4.71238898038469 6.28318530717959],...
        'XTickLabel',{'0','$\frac{\pi}{2}$','$\pi$','$\frac{3\pi}{2}$','$2 \pi$'},...
        'YMinorGrid','on','ZGrid','on');
    % Create legend
    legend1 = legend(axes1,'show');
    set(legend1,...
        'Position',[0.0828722015781729 0.82523809387571 0.339263807024274 0.167380955105736],...
        'Interpreter','latex',...
        'EdgeColor','none');
end