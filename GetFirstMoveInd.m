function firstMoveInd = GetFirstMoveInd(va)
    % GETFIRSTMOVEIND takes a data array of velocity*accel which can be
    % interpretted as intent of movement.  It finds the last point in time
    % where va < 0.005  before va is maximized
    

    [maxval, maxind] = max(va(1:end-2000));
    inds = va(1:maxind) < 0.005;
    firstMoveInd = find(inds, 1, 'last');
    
end