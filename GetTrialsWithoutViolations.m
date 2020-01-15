function trials_wov = GetTrialsWithoutViolations(trials, targetDirNum, dampNum)
    
    trials_wov = [];
    
    % Get subset of trials that match with inputted target direction number
    % and damp number, and had group number >= 4
    groupNum_all = [trials.GroupNumber];    
    targetNum_all = [trials.TargetDirNum];
    dampNum_all = [trials.DampingNumber];
    
    groupNum_inds = (groupNum_all >= 4);
    targetNum_inds = (targetNum_all == targetDirNum);
    dampNum_inds = (dampNum_all == dampNum);
    
    set_inds = (groupNum_inds & targetNum_inds & dampNum_inds);
    trials_set = trials(set_inds);
    
    % Get x data
    nSamplesSet = [trials_set.nSamples];
    nSamples = min(nSamplesSet);
    
    nTrials = length(trials_set);
    x = zeros(nSamples, nTrials);
    
    for i = 1:nTrials
        x(:,i) = trials_set(i).Data.EndEffPos_FromJA(1:nSamples);
    end
    
    x = GetDataSetStats(x);
    
    % Test each trial for violation of +/- 3 std bound
    for i = 1:nTrials
        x_temp = trials_set(i).Data.EndEffPos_FromJA(1:nSamples);
        
        ub_violations = sum(x_temp > x.ub_);
        lb_violations = sum(x_temp < x.lb_);
        total_violations = ub_violations + lb_violations;
        
        if (total_violations == 0)
            trials_wov = [trials_wov; trials_set(i)];
        end
    end
    
end