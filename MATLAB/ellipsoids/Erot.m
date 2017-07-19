function result = Erot(w,I)
% Returns the translationnal energy of a particle of velocity v and mass
% m.

    w = transpose(w); % rotation vector

    I_ = eye(3); % moment of inertia tensor
    for coord = 1:3
        I_(coord,coord) = I(coord);
    end

    result = 1/2 * w * I_ * transpose(w);
end