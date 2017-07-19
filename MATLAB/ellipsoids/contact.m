function result = contact(E1,r1,q1,g1,E2,r2,q2,g2)
% Returns the contact points of two overlapping ellipsoids of belonging
% matrix E1 and E2 and ever-belonging balls' radius g1 and g2, whose
% centers are r1 and r2 and orientations are described by q1 and q2.
% These contact points are determined with a rescaling method inspired by
% part III.A of Donev et al., Phys. Rev. E 2007 based on Perram and
% Wertheim, J. Comput. Phys. 1985.

    epsilon_dist = 1e-4; % precision on distances

    c = distellipsoid(E1,r1,q1,g1,E2,r2,q2,g2); % closest points between the ellipsoids
    
    if c ~= -1 % the ellipsoids do not overlap or are almost touching
        result = -1; % the resulting potential will be zero, we are not considering any interaction
        return
    end
    
    scaling_min = 1/2; % minimum rescaling searched
    scaling = [scaling_min 1]; % interval of scaling searched
    
    dist = -1; % distance between rescaled ellipsoids
    
    while dist > epsilon_dist || dist < 0 % rescaled ellipsoids are too far from each other
        
        scaling_ = sum(scaling)/2; % rescaling factor
        c = distellipsoid(dilat(1/scaling_)*E1,r1,q1,g1,dilat(1/scaling_)*E2,r2,q2,g2); % closest points between the rescaled ellipsoids
        if c == -1 % the rescaling factor is too high, rescaled ellipsoids still overlap
            scaling(2) = scaling_;
        else % the rescaling factor is too low, rescaled ellipsoids are too far from each other
            dist = norm(c(1,:)-c(2,:)); % distance between the rescaled ellipsoids
            scaling(1) = scaling_;
        end
    end
    
    scaling = sum(scaling)/2; % computed scaling for the rescaled ellipsoids to be extrernally tangent
    
    c(1,:) = transpose(quat2vec(trans(r1)*quat2rot(q1)*dilat(1/sqrt(scaling))*transpose(quat2rot(q1))*trans(-r1)*transpose(vec2hom(c(1,:))))); % corresponding contact point on ellipsoid 1
    c(2,:) = transpose(quat2vec(trans(r2)*quat2rot(q2)*dilat(1/sqrt(scaling))*transpose(quat2rot(q2))*trans(-r2)*transpose(vec2hom(c(2,:))))); % corresponding contact point on ellipsoid 2

    result = c;

end