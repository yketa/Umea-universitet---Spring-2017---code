function result = fdis(kd,contacts,r1,v1,q1,qp1,r2,v2,q2,qp2)
% Returns the dissipation force exerted on particle 1 by particle 2.

    v1_ = v1 + cross(w(q1,qp1),contacts(1,:)-r1); % velocity of the contact point of particle 1
    v2_ = v2 + cross(w(q2,qp2),contacts(2,:)-r2); % velocity of the contact point of particle 2
    
    result = - kd * (v1_ - v2_);
    
end