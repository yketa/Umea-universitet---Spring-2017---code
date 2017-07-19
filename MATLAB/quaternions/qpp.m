function result = qpp(q,qp,M,I)
% Returns the second derivative of quaternion q, according to q, its first
% derivative qp, the moment M exerted on the particle and the particle's
% moment of inertia.

    result = qprod(q,vec2quat(M/(2*I)) - scal2quat(qnorm(qp)));
end