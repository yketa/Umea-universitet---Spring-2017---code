function result = Erot(w,I)
% Returns the translationnal energy of a particle of velocity v and mass
% m.

    result = 1/2 * I * sum(w.^2,2);
end