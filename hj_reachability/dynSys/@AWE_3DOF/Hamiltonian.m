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