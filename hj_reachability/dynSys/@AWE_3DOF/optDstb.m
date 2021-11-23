function dOpt = optDstb(obj, ~, y, deriv, dMode, ~)
% dOpt = optDstb(obj, t, y, deriv, ~, ~)

if nargin < 5
  dMode = 'max';
end

if ~(strcmp(dMode, 'max') || strcmp(dMode, 'min'))
  error('dMode must be ''max'' or ''min''!')
end


if ~iscell(deriv)
  deriv = num2cell(deriv);
end
dOpt = -sign(deriv{4}) .* obj.max_tether_diff_dot;

end