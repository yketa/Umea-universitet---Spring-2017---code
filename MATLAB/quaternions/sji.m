function result = sji(ri,di,rj,dj)
% Returns the moment arm point from the center of particle i to point of
% contact with particle j.

    dij = (di + dj)/2; % average diameter
    Ri = di/2; % radius of particle i

    result = (rj - ri)*(Ri/dij);
    
end