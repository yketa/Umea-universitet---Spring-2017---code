function result = qunit(q)
% Returns the unit quaternion associated to q.

    if qnorm(q) ~= 0 % the quaternion is not 0
        result = q/sqrt(qnorm(q));
    else % the quaternion is 0
        result = q;
        
end