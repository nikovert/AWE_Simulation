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

function [output] = add_cellMatrix(mat1,mat2)
%MULT_CELLMATRIX Multiples two Matrices
[r1, c1] = size(mat1);
[r2, c2] = size(mat2);
% if nargin < 3  
%     r3 = r2;
%     c3 = c2;
%     mat3 = 0;
% else
%     [r3, c3] = size(mat3);
% end
% assert((r1 == r2 == r3 && c1 == c2 == c3) || (r1 == 1 && c1 == 1) || (r2 == 1 && c2 == 1))
assert(r1 == r2 && c1 == c2)
if ~iscell(mat1)
    mat1 = num2cell(mat1, 2)';
end

if ~iscell(mat2)
    mat2 = num2cell(mat2, 2)';
end

output = cell(r1, c1);
for column = 1:c2
    for row = 1:r1
        output{row, column} = mat1{row, column} + mat2{row, column};
    end
end
end

