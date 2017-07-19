function result = fdis(kd,ri,vi,wi,di,rj,vj,wj,dj)
% Returns the dissipation force exerted on particle i by particle j.

    rij = ri - rj;
    Rij = norm(rij);
    dij = (di+dj)/2;
    
    if Rij < dij % particles touch
        result = -kd*(vi - vj + cross(wi,sji(ri,di,rj,dj)) - cross(wj,sji(rj,dj,ri,di)));
    else % particles do not touch
        result = [0 0 0];
    end
    
end