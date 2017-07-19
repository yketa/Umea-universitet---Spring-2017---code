function result = distellipsoid(E1,r1,q1,g1,E2,r2,q2,g2)
% Returns the closest points between 2 ellispoids of belonging matrix E1
% and E2 and ever-belonging balls' radius g1 and g2, whose centers are r1
% and r2 and orientations are described by q1 and q2.
% This algorithm relies on the method described by Lin and Han, SIAM
% Journal of Optimisation, 2002.

    epsilon_dist = 1e-4; % precision on distances
    epsilon_ang = 1e-3; % precision on angles
    
    % determination of the intersection of r1r2 with both ellipsoids
    
    t1 = [0 1]; % interval of search of the intersection with ellipsoid 1
    t2 = [0 1]; % interval of search of the intersection with ellipsoid 2
    
    while norm((1-t1(1))*r1 + t1(1)*r2 - (1-t1(2))*r1 - t1(2)*r2) > epsilon_dist && norm((1-t2(1))*r1 + t2(1)*r2 - (1-t2(2))*r1 - t2(2)*r2) > epsilon_dist % the computed intersections are not close enough to the ellipsoids
        
        t1_ = sum(t1)/2;
        if inellipsoid(E1,r1,q1,(1-t1(1))*r1 + t1(1)*r2) + inellipsoid(E1,r1,q1,(1-t1_)*r1 + t1_*r2) == 2 % the interval [t1(1);t1_] is in the ellipsoid 1
            t1(1) = t1_;
        else
            t1(2) = t1_;
        end
        
        t2_ = sum(t2)/2;
        if inellipsoid(E2,r2,q2,(1-t2(1))*r1 + t2(1)*r2) + inellipsoid(E2,r2,q2,(1-t2_)*r1 + t2_*r2) == 0 % the interval [t2(1);t2_] is not in the ellipsoid 2
            t2(1) = t2_;
        else
            t2(2) = t2_;
        end
        
    end
    
    t1 = sum(t1)/2; % computed intersection with ellipsoid 1
    t2 = sum(t2)/2; % computed intersection with ellipsoid 2
    
    c1 = (1-t1)*r1 + t1*r2; % corresponding point on ellipsoid 1
    c2 = (1-t2)*r1 + t2*r2; % corresponding point on ellipsoid 2
    
    n1 = quat2vec(belonging(E1,r1,q1)*transpose(vec2hom(c1))); % outward facing orthogonal vector to the surface of ellipsoid 1 at c1
    n2 = quat2vec(belonging(E2,r2,q2)*transpose(vec2hom(c2))); % outward facing orthogonal vector to the surface of ellipsoid 2 at c2
    
    if sum(c1-c2) < epsilon_dist % the points are almost touching
        result = [c1;c2];
        return
    elseif t2 <= t1 % the ellipsoids are overlapping
        result = -1;
        return
    elseif angle(c2 - c1,n1) < epsilon_ang && angle(c1 - c2,n2) < epsilon_ang % the line between the computed points is aligned with the surface vectors, then the minimum distance has been found
        result = [c1;c2];
        return
    end
    
    % determination of the distance between ellipsoids
    
    while t2 > t1 && (angle(c2 - c1,n1) > epsilon_ang || angle(c1 - c2,n2) > epsilon_ang)
        
        r1_ = c1 - g1*transpose(n1); % center of the ball constructed in ellipsoid 1
        r2_ = c2 - g2*transpose(n2); % center of the ball constructed in ellipsoid 2
        
        t1 = [0 1]; % interval of search of the intersection with ellipsoid 1
        t2 = [0 1]; % interval of search of the intersection with ellipsoid 2

        while norm((1-t1(1))*r1_ + t1(1)*r2_ - (1-t1(2))*r1_ - t1(2)*r2_) > epsilon_dist && norm((1-t2(1))*r1_ + t2(1)*r2_ - (1-t2(2))*r1_ - t2(2)*r2_) > epsilon_dist % the computed intersections are not close enough to the ellipsoids

            t1_ = sum(t1)/2;
            if inellipsoid(E1,r1,q1,(1-t1(1))*r1 + t1(1)*r2) + inellipsoid(E1,r1,q1,(1-t1_)*r1 + t1_*r2) == 2  % the interval [t1(1);t1_] is in the ellipsoid 1
                t1(1) = t1_;
            else
                t1(2) = t1_;
            end

            t2_ = sum(t2)/2;
            if inellipsoid(E2,r2,q2,(1-t2(1))*r1 + t2(1)*r2) + inellipsoid(E2,r2,q2,(1-t2_)*r1 + t2_*r2) == 2 % the interval [t2(1);t2_] is not in the ellipsoid 2
                t2(2) = t2_;
            else
                t2(1) = t2_;
            end

        end

        t1 = sum(t1)/2; % computed intersection with ellipsoid 1
        t2 = sum(t2)/2; % computed intersection with ellipsoid 2
        
        c1 = (1-t1)*r1_ + t1*r2_; % corresponding point on ellipsoid 1
        c2 = (1-t2)*r1_ + t2*r2_; % corresponding point on ellipsoid 2
        
        n1 = quat2vec(belonging(E1,r1,q1)*transpose(vec2hom(c1))); % outward facing orthogonal vector to the surface of ellipsoid 1 at c1
        n2 = quat2vec(belonging(E2,r2,q2)*transpose(vec2hom(c2))); % outward facing orthogonal vector to the surface of ellipsoid 2 at c2
        
    end
    
    if sum(c1-c2) < epsilon_dist  % the points are almost touching
        result = [c1;c2];
        return
    elseif t2 <= t1 % the ellipsoids are overlapping
        result = -1;
        return
    elseif angle(c2 - c1,n1) < epsilon_ang && angle(c1 - c2,n2) < epsilon_ang % the line between the computed points is aligned with the surface vectors, then the minimum distance has been found
        result = [c1;c2];
        return
    end
    
end