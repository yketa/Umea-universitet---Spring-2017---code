function result = qrot(rot)
% Returns the quaternion associated to the rotations rot(1) to rot(end)
% in this order.
% rot(i) have to be submitted as [ax,ang] where ax = 1, 2 or 3
% corresponding to the axis x, y or z, and ang corresponding to the angle
% of the rotation along this axis.

    result = [0 0 0 1]; % quaternion
    
    for rotation = 1:length(rot)
        
        axis = [0 0 0]; % axis of the rotation
        axis(rot(rotation,1)) = 1;
        
        angle = rot(rotation,end); % angle of the rotation
        
        result = qprod(vecscal2quat(axis*sin(angle/2),cos(angle/2)),result);
        
    end