function res = FilterAndAnalyzeData(trials, TargetDirNum, DampNum, varargin)
    res = struct;
    % Parse arguments
    p = inputParser;
    defaultVal = -1;
    addParameter(p, 'SubjectNum', defaultVal);
    
    parse(p, varargin{:});
    
    % Get inputted values to filter by
    subjectNumInput = p.Results.SubjectNum;
    targetDirNumInput = TargetDirNum;
    dampNumInput = DampNum;
    
    % DefaultInds
    groupNumInds = ([trials.GroupNumber] >= 4); % default
    subjectNumInds = logical(ones(1, length(trials)));
    targetDirNumInds = logical(ones(1, length(trials)));
    dampNumInds = logical(ones(1, length(trials)));
    nSamplesInds = ([trials.nSamples] > 2000); % all trials should be at least 2 sec long, thus at least 2000 samples in length
    
    % Refine Inds To Filter With Inputs
    if (subjectNumInput ~= defaultVal)
        subjectNumInds = ([trials.SubjectNumber] == subjectNumInput);
    end
    
    if (targetDirNumInput ~= defaultVal)
        targetDirNumInds = ([trials.TargetDirNum] == targetDirNumInput);
    end
    
    if (dampNumInput ~= defaultVal)
        dampNumInds = ([trials.DampingNumber] == dampNumInput);
    end
    
    % Combine Inds and Filter
    filterInds = (groupNumInds & subjectNumInds & targetDirNumInds & dampNumInds & nSamplesInds);
    trials_set = trials(filterInds);
    
    % Exit function if no trials match these criteria
    if (numel(trials_set) == 0)
        fprintf('No trials with group number >=4, target direction: %d, damping number: %d.\n', ...
                targetDirNum, ...
                dampNum);
        res = [];
        return;
    end
    
    % Find min nSamples in set
    nSamplesSet = [trials_set.nSamples];
    nSamples = min(nSamplesSet);
    res.nSamples = nSamples;
    
    % Get all end effector positions calc'ed from JA, get array of trials
    % without violations (wov) and with violations (wv) of the +/- 3 std of
    % regular position movement
    nTrials = length(trials_set);
    x = zeros(nSamples, nTrials);
    for i = 1:length(trials_set)
        x(:,i) = trials_set(i).Data.EndEffPos_FromJA(1:nSamples);
    end
    
    x_mean = mean(x, 2);
    x_std = std(x, [], 2);
    x_ub = x_mean + 3*x_std;
    x_lb = x_mean - 3*x_std;
    
    inds_wov = logical(zeros(nTrials, 1));
    inds_wv = logical(zeros(nTrials, 1));
    for i = 1:nTrials
        x_col = x(:,i);
        ub_viol = x_col > x_ub;
        lb_viol = x_col < x_lb;
        
        if ((sum(ub_viol) > 0) || (sum(lb_viol) > 0))   % violation of ub or lb occured
            inds_wv(i) = true;
        else                                            % no violatons occured
            inds_wov(i) = true;        
        end
    end
    
    trials_wov = trials_set(inds_wov);
    trials_wv = trials_set(inds_wv);
    res.wov.trials = trials_wov;
    res.wv.trials = trials_wv;

    % Get data sets
    nTrialsWov = length(trials_wov);
    nTrialsWv = length(trials_wv);
    
    x_wov   = zeros(nSamples, nTrialsWov);
    x_wv    = zeros(nSamples, nTrialsWv);
    va_wov  = zeros(nSamples, nTrialsWov);
    va_wv   = zeros(nSamples, nTrialsWv);
    damp_wov    = zeros(nSamples, nTrialsWov);
    damp_wv     = zeros(nSamples, nTrialsWv);
    xdot_wov    = zeros(nSamples, nTrialsWov);
    xdot_wv     = zeros(nSamples, nTrialsWv);
    force_wov   = zeros(nSamples, nTrialsWov);
    force_wv    = zeros(nSamples, nTrialsWv);
    
    for i = 1:nTrialsWov
        x_wov(:,i)      = trials_wov(i).Data.EndEffPos_FromJA(1:nSamples);
        va_wov(:,i)     = trials_wov(i).Data.va(1:nSamples);
        damp_wov(:,i)   = trials_wov(i).Data.Damping(1:nSamples);
        xdot_wov(:,i)   = trials_wov(i).Data.xdot(1:nSamples);
        force_wov(:,i)      = trials_wov(i).Data.Force(1:nSamples, 1);
    end
    
    res.wov.x = GetDataSetStats(x_wov);
    res.wov.va = GetDataSetStats(va_wov);
    res.wov.damp = GetDataSetStats(damp_wov);
    res.wov.xdot = GetDataSetStats(xdot_wov);
    res.wov.force = GetDataSetStats(force_wov);
    
    for i = 1:nTrialsWv
        x_wv(:,i)      = trials_wv(i).Data.EndEffPos_FromJA(1:nSamples);
        va_wv(:,i)     = trials_wv(i).Data.va(1:nSamples);
        damp_wv(:,i)   = trials_wv(i).Data.Damping(1:nSamples);
        xdot_wv(:,i)   = trials_wv(i).Data.xdot(1:nSamples);
        force_wv(:,i)      = trials_wv(i).Data.Force(1:nSamples, 1);
    end
    
    res.wv.x = GetDataSetStats(x_wv);
    res.wv.va = GetDataSetStats(va_wv);
    res.wv.damp = GetDataSetStats(damp_wv);
    res.wv.xdot = GetDataSetStats(xdot_wv);
    res.wv.force = GetDataSetStats(force_wv);
    

    % Find when movement starts -- wov
    firstMoveInds = zeros(nTrialsWov,1);
    restOfArray = zeros(nTrialsWov, 1);
    for i = 1:nTrialsWov
        va_temp = trials_wov(i).Data.va;
        firstMoveInds(i) = GetFirstMoveInd(va_temp);
        restOfArray(i) = length(va_temp) - firstMoveInds(i);
    end
    
    minRestOfArray = min(restOfArray);
    x_mm = [];          % movement matching
    va_mm = [];
    damp_mm = [];
    force_mm = [];
    for i = 1:nTrialsWov
        firstInd = firstMoveInds(i);% - minFirstMoveInd + 1;
        lastInd = firstMoveInds(i) + minRestOfArray;
        x_mm = [x_mm, trials_wov(i).Data.EndEffPos_FromJA(firstInd:lastInd)];
        va_mm = [va_mm, trials_wov(i).Data.va(firstInd:lastInd)];
        damp_mm = [damp_mm, trials_wov(i).Data.Damping(firstInd:lastInd)];
        force_mm = [force_mm, trials_wov(i).Data.Force(firstInd:lastInd, 1)];
    end
    res.wov.x_mm = GetDataSetStats(x_mm);
    res.wov.va_mm = GetDataSetStats(va_mm);
    res.wov.damp_mm = GetDataSetStats(damp_mm);
    res.wov.force_mm = GetDataSetStats(force_mm);
   
    
    

    % Record where the data that was analyzed came from
    res.SubjectNumbers = unique([trials_set.SubjectNumber]);
    res.TargetDirNum = TargetDirNum;
    res.DampingNum = DampNum;
    res.TargetDirText = trials_set(1).TargetDirText;
    res.DampingText = trials_set(1).DampingText;
    
    % Calculate Overshoot
    res.PercentOvershoot = GetPercentOvershootStats(trials_wov, nSamples);
    
    % Calculate max velocity
    res.XdotMax = GetMaxVelocityStats(trials_wov, nSamples);
    
end