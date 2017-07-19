function result = qnorm(quat)
% Returns the norm a quaternion.

    result = sum(quat.^2,2);
end