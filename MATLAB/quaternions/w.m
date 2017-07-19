function result = w(q,qp)
% Returns the angular velocity associated to a quatrernion and its
% first derivative.

    result = 2*quat2vec(qprod(qconj(q),qp));
end