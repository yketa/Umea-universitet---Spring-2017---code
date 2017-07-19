function result = trans(vec)
% Returns the translation matrix of vector vec acting on homogeneous
% coordinates.

    result = eye(4);
    result(1:3,4) = transpose(vec);
    
end