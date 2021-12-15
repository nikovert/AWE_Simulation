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

function x = checkInterpInput(g, x)
if size(x, 2) ~= g.dim
  if size(x, 1) == g.dim
    % Take transpose if number of input rows is same as grid dimension
    x = x';
  else
    error('Input points must have the same dimension as grid!')
  end
end
end