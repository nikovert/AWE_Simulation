function [output] = cross_cellVec(vec1,vec2)
%cross_cellVec  Vector cross product.
%   C = CROSS(A,B) returns the cross product of the vectors
%   A and B.  That is, C = A x B.  A and B must be 3 element
%   vectors.
%
%   C = CROSS(A,B) returns the cross product of A and B along the
%   first dimension of length 3.
    
    % Calculate cross product
    output = {vec1{2}.*vec2{3}-vec1{3}.*vec2{2}; ...
                vec1{3}.*vec2{1}-vec1{1}.*vec2{3};...
                vec1{1}.*vec2{2}-vec1{2}.*vec2{1}};
end

