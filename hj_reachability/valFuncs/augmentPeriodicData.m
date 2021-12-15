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

function [g, data] = augmentPeriodicData(g, data)

%% Dealing with periodicity
for i = 1:g.dim
  if isfield(g, 'bdry') && isequal(g.bdry{i}, @addGhostPeriodic)
    % Grid points
    g.vs{i} = cat(1, g.vs{i}, g.vs{i}(end) + g.dx(i));
    
    % Input data; eg. data = cat(:, data, data(:,:,1))
    colons1 = repmat({':'}, 1, g.dim);
    colons1{i} = 1;
    cat_argin = {i; data; data(colons1{:})};
    data = cat(cat_argin{:});
  end
end

end