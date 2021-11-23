function [output] = transpose_cellMatrix(input)
%TRANSPOSE_CELL transposes a vecotr contained in a cell
[r, c] = size(input);
output = cell(c, r);
for row = 1:r
    for column = 1:c
        output{column, row} = input{row, column};
    end
end
end
