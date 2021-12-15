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



%% Add any neccessary toolboxes to the matlab path
currentFile = mfilename( 'fullpath' );
[pathstr,~,~] = fileparts( currentFile );
addpath(genpath([pathstr, '/../../Add-Ons/level-set-methods-toolbox/Kernel/']));
addpath([pathstr, '/visualization'])
addpath([pathstr, '/valFuncs'])
addpath([pathstr, '/grids'])
addpath([pathstr, '/dynSys'])