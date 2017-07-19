function result = Etrans(v,m)
% Returns the translationnal energy of a particle of velocity v and mass
% m.

    result = 1/2 * m * sum(v.^2,2);
end