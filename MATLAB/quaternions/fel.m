function result = fel(ke,ri,di,rj,dj)
% Returns the elastic force exerted on particle i by particle j.

    rij = ri - rj;
    Rij = norm(rij);
    dij = (di+dj)/2;
    
    if Rij < dij % particles touch
        result = ke/(dij*Rij) * (1 - Rij/dij) * rij;
    else % particles do not touch
        result = [0 0 0];
    end
    
end