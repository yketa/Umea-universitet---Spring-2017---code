function result = dilat(a)
% Returns dilatation matrix of factor a acting on homogeneous coordinates.

    result = a*eye(4);
    result(4,4) = 1;
end