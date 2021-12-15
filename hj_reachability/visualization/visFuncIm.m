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
% :Author: Mo Chen (mochen@sfu.ca)

function [h]= visFuncIm(gPlot,dataPlot,color,alpha,fignr)

if nargin < 5
    fignr = 1;
end
figure(fignr);

if gPlot.dim<2
    h = plot(gPlot.xs{1}, squeeze(dataPlot),...
        'LineWidth',2);
    h.Color = color;
elseif gPlot.dim==2
    h = surf(gPlot.xs{1}, gPlot.xs{2}, dataPlot);
    h.EdgeColor = 'none';
    h.FaceColor = color;
    h.FaceAlpha = alpha;
    h.FaceLighting = 'phong';
else
    error('Can not plot in more than 3D!')
end

end
