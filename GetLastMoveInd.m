function lastMoveInd = GetLastMoveInd(va)

    % GETLASTMOVEIND takes a data array of velocity*accel which can be
    % interpretted as intent of movement.  It finds the last point in time
    % where abs(va > 0.001)  before va is maximized
    
    inds = abs(va) > 0.001;
    lastMoveInd = find(inds, 1, 'last');
    

end