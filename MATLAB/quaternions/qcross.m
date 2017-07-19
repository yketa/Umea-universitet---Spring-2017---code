function result = qcross(quat1,quat2)
% Returns the cross product of vectors or quaternions quat1 and quat2.

    result = cross(quat1(1:3),quat2(1:3));
end