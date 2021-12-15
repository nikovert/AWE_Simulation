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

function hamValue = Hamiltonian(obj, t, ~, deriv, schemeData)
% function hamValue = genericHam(t, data, deriv, schemeData)

%% Input unpacking
if ~isfield(schemeData, 'uMode')
    schemeData.uMode = 'min';
end

if ~isfield(schemeData, 'dMode')
    schemeData.dMode = 'max';
end

if ~isfield(schemeData, 'tMode')
    schemeData.tMode = 'backward';
end

%% Plug optimal control into dynamics to compute Hamiltonian
obj.optCtrlMode = 'normal';
uOpt = obj.optCtrl(t, schemeData.grid.xs, deriv, schemeData.uMode);
dOpt = obj.optDstb(t, schemeData.grid.xs, deriv, schemeData.dMode);
dx = obj.dynamics(t, schemeData.grid.xs, uOpt, dOpt);
hamValue = 0;
for i = 1:obj.nx
  hamValue = hamValue + deriv{i}.*dx{i};
end

%% Negate hamValue if backward reachable set
if strcmp(schemeData.tMode, 'backward')
    hamValue = -hamValue;
end

if isfield(schemeData, 'side')
    if strcmp(schemeData.side, 'upper')
        hamValue = -hamValue;
    end
end

end