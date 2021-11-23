function [output] = norm_cellVec(input)
    output = sqrt(input{1}.^2 + input{2}.^2 + input{3}.^2);
end



