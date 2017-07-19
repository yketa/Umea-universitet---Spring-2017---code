function result = qpp(q,qp,M,I)
% Returns the second derivative of quaternion q, according to q, its first
% derivative qp, the moment M exerted on the particle and the particle's
% moments of inertia.

    omega = transpose(w(q,qp)); % rotation vector
    M = transpose(M); % sum of moments
    
    I_ = eye(3); % moment of inertia tensor
    I_inv = eye(3); % inverse of the moment of inertia tensor
    for coord = 1:3
        I_(coord,coord) = I(coord);
        I_inv(coord,coord) = 1/I(coord);
    end

    result = qprod(q,vec2quat(transpose(1/2 * I_inv * (M - cross(omega,I_*omega)))) - scal2quat(qnorm(qp)));
    
end