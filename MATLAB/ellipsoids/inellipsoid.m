function result = inellipsoid(E,r,q,vec)
% Returns 1 if r in is in the ellipsoid of belonging matrix E whose center
% is r and orientation is described by q, 0 otherwise.

    result = vec2hom(vec)*belonging(E,r,q)*transpose(vec2hom(vec)) <= 0;
    
end