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

function alpha = partialFunc(obj, t, data, derivMin, derivMax, schemeData, dim)
% alpha = genericPartial(t, data, derivMin, derivMax, schemeData, dim)

g = schemeData.grid;
dynSys = schemeData.dynSys;

if ~isfield(schemeData, 'uMode')
  schemeData.uMode = 'min';
end

if ~isfield(schemeData, 'dMode')
  schemeData.dMode = 'min';
end

%% Compute control
% Optimal control assuming maximum deriv
if dim > 3 && dim < 7
    dynSys.optCtrlMode = 'derivMax';
    uU = dynSys.optCtrl(t, g.xs, derivMax, schemeData.uMode);

    % Optimal control assuming minimum deriv
    dynSys.optCtrlMode = 'derivMin';
    uL = dynSys.optCtrl(t, g.xs, derivMin, schemeData.uMode);

    %% Compute disturbance
    dU = {0; 0; 0; 0};
    dL = {0; 0; 0; 0};

    %% Compute alpha
    dxUU = dynSys.dynamics(t, schemeData.grid.xs, uU, dU);
    dxLL = dynSys.dynamics(t, schemeData.grid.xs, uL, dL);
    alpha = max(abs(dxUU{dim}), abs(dxLL{dim}));
elseif dim == 7 || dim < 4
    uU = {0;0};
    uL = {0;0};
    dU = obj.optDstb(t, schemeData.grid.xs, derivMax, schemeData.dMode);
    dL = obj.optDstb(t, schemeData.grid.xs, derivMin, schemeData.dMode);
    %% Compute alpha
    dxUU = dynSys.dynamics(t, schemeData.grid.xs, uU, dU);
    dxLL = dynSys.dynamics(t, schemeData.grid.xs, uL, dL);
    alpha = max(abs(dxUU{dim}), abs(dxLL{dim}));
else
    if isempty(obj.long_dot)
        dynSys.dynamics(t, schemeData.grid.xs, {0,0}, 0);
    end
    switch dim
        case 1
            alpha = abs(obj.long_dot);
        case 2
            alpha = abs(obj.lat_dot);
        case 3
            alpha = abs(obj.h_tau_dot);
        case 7
            n_t_p = obj.T.np; % amount of particles
            l_s_dot = obj.v_ro_set/(n_t_p+1);
            alpha = abs(l_s_dot).*ones(schemeData.grid.shape);
        otherwise 
            error('Ups')
    end
end
