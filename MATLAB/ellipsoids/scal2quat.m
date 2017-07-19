function result = scal2quat(scal)
% Returns the quaternion associated to a scalar.

    result = cat(2,[0 0 0],scal);
end