function result = quat2rot(q)
% Returns the rotation matrix associated to the quaternion q acting on
% homogeneous coordinates.

    cross = [0 -q(3) q(2);q(3) 0 -q(1);-q(2) q(1) 0]; % matrix associated to the cross product with vector part of q
    mat = (cross + q(4)*eye(3))^2 + transpose(quat2vec(q))*quat2vec(q); % rotation matrix
    
    result = eye(4);
    result(1:3,1:3) = mat;
    
end