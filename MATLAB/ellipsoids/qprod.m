function result = qprod(quat1,quat2)
% Returns the products of quaternions
% quat1 = [vec1,scal1] and quat2 = [vec2,scal2]
% such as quat = scal + i vec.1 + j vec.2 + k vec.3

    vector1 = cross(quat1(1:3),quat2(1:3)); % first term in the sum of vectors
    vector2 = quat1(4)*quat2(1:3); % second term in the sum of vectors
    vector3 = quat2(4)*quat1(1:3); % third term in the sum vectors
    scalar = quat1(4)*quat2(4) - dot(quat1(1:3),quat2(1:3)); % scalar part of the product

    result = cat(2,sum([vector1;vector2;vector3],1),scalar); % product of quaternions
end