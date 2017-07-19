function result = belonging(E,r,q)
% Returns the adjusted belonging matrix of ellipsoid of belonging matrix E
% whose center is r and orientation is described by q.

    translation = trans(-r); % translation matrix
    rotation = quat2rot(q); % rotation matrix

    result = transpose(translation)*rotation*E*transpose(rotation)*translation; % adjusted belonging matrix
    
end