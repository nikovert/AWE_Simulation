function [output] = element_div_cellMatrix(mat1,mat2)
%MULT_CELLMATRIX Multiples two Matrices

if ~iscell(mat2)
    mat2 = {mat2};
end

[r1, c1] = size(mat1);
[r2, c2] = size(mat2);

assert(r1 == r2 && c1 == c2 || (r1 == 1 && c1 == 1) || (r2 == 1 && c2 == 1))
if (r1 == 1 && c1 == 1)
    output = cell(r2, c2);
    for column = 1:c2
        for row = 1:r2
            output{row, column} = mat1{1, 1} ./ mat2{row, column};
        end
    end
elseif (r2 == 1 && c2 == 1)
    output = cell(r1, c1);
    for column = 1:c1
        for row = 1:r1
            output{row, column} = mat1{row, column} ./ mat2{1, 1};
        end
    end
else
    output = cell(r1, c2);
    for column = 1:c2
        for row = 1:r1
            output{row, column} = mat1{row, column} ./ mat2{row, column};
        end
    end
end
end

