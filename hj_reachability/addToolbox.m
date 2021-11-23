%% Add any neccessary toolboxes to the matlab path
currentFile = mfilename( 'fullpath' );
[pathstr,~,~] = fileparts( currentFile );
addpath(genpath([pathstr, '/../../Add-Ons/level-set-methods-toolbox/Kernel/']));
addpath([pathstr, '/visualization'])
addpath([pathstr, '/valFuncs'])
addpath([pathstr, '/grids'])
addpath([pathstr, '/dynSys'])