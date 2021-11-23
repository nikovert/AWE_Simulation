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
    dU = 0;
    dL = 0;

    %% Compute alpha
    dxUU = dynSys.dynamics(t, schemeData.grid.xs, uU, dU);
    dxLL = dynSys.dynamics(t, schemeData.grid.xs, uL, dL);
    alpha = max(abs(dxUU{dim}), abs(dxLL{dim}));
elseif dim == 7
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
