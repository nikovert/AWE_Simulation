function [output] = mult_cellMatrix(mat1,mat2)
%MULT_CELLMATRIX Multiples two Matrices
[r1, c1] = size(mat1);
[r2, c2] = size(mat2);
if ~iscell(mat1)
    mat1 = num2cell(mat1);
end
if ~iscell(mat2)
    mat2 = num2cell(mat2);
end
assert((r2 == c1) || (r1 == 1 && c1 == 1) || (r2 == 1 && c2 == 1))

if (r1 == 1 && c1 == 1)
    output = cell(r2, c2);
    for column = 1:c2
        for row = 1:r2
            output{row, column} = mat2{row, column} .* mat1{1, 1};
        end
    end
elseif (r2 == 1 && c2 == 1)
    output = cell(r1, c1);
    for column = 1:c1
        for row = 1:r1
            output{row, column} = mat1{row, column} .* mat2{1, 1};
        end
    end
else
    output = cell(r1, c2);
    for column = 1:c2
        for row = 1:r1
            val = 0;
            for index = 1:r2
                val = val + mat1{row, index} .* mat2{index, column};
            end
            output{row, column} = val;
        end
    end
end
end

