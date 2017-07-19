function result = fel(ke,d,overlap,contact_point,E,r,q)
% Returns the elastic force exerted at contact_point on a particle of
% belonging matrix E whose center is r and orientation is described by q.

    direction = quat2vec(belonging(E,r,q)*transpose(vec2hom(contact_point))); % direction of the force
    direction = direction/norm(direction); % the direction has to be an unit vector

    result = - ke/d * (1 - (d -overlap)/d) * transpose(direction);
    
end