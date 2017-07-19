function result = angle(a,b)
% Returns the angle between vectors a and b.

    result = atan2(norm(cross(a,b)),dot(a,b));
end