function result = Mdis(kd,ri,vi,wi,di,rj,vj,wj,dj)
% Returns the moment of the dissipation force exerted on the particle i
% on the particle j.

    result = cross(sji(ri,di,rj,dj),fdis(kd,ri,vi,wi,di,rj,vj,wj,dj));
end