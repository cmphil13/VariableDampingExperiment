function res = GetDataSetStats(data)
    % GETDATASETSTATS takes a data set of nSamples x nTrials size, and returns
    % the group's mean, std, 3 std upper bound (ub), 3 std lower bound (lb)
    % in a struct

    res = struct;
    
    res.data_ = data;
    res.mean_ = mean(data, 2);
    res.std_ = std(data, [], 2);
    res.ub_ = res.mean_ + 3*res.std_;
    res.lb_ = res.mean_ - 3*res.std_;
    res.nSamples_ = size(data, 1);
    res.nTrials_ = size(data, 2);

end