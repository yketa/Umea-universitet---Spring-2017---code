function result = qconj(quat)
% Returns the conjugate of a quaternion.

    result = cat(2,-quat(1:3),quat(4));
end