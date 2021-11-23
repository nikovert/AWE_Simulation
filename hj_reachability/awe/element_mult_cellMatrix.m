function [output] = element_mult_cellMatrix(mat1,mat2)
%MULT_CELLMATRIX Multiples two Matrices
[r1, c1] = size(mat1);
[r2, c2] = size(mat2);

assert(r1 == r2 && c1 == c2)
output = cell(r1, c2);
for column = 1:c2
    for row = 1:r1
        output{row, column} = mat1{row, column} .* mat2{row, column};
    end
end
end

