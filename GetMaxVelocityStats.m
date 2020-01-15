function res = GetMaxVelocityStats(trials_wov, nSamples)

    % GETMAXVELOCITYSTATS takes in a the whole trial struct of trials
    % without violations and nSamples.  It returns a struct array of the
    % individual trial max velocity, max, mean, std for
    % each subject.
    
    res = [];
    
    subjectNumbersAll = [trials_wov.SubjectNumber];
    subjectNumbers = unique(subjectNumbersAll);
    
    for i = 1:length(subjectNumbers)
        subjectNumber = subjectNumbers(i);
        subjectInds = subjectNumbersAll == subjectNumber;
        
        trials_subject = trials_wov(subjectInds);
        nTrialsSubject = length(trials_subject);
        xdot = zeros(nSamples, nTrialsSubject);
        for j = 1:nTrialsSubject
            xdot(:,j) = abs(trials_subject(j).Data.xdot(1:nSamples));
        end
        
        xdot_max_struct = struct;       
        
        xdot_max_struct.max_xdot = max(xdot, [], 1);
        xdot_max_struct.max_xdot_mean = mean(xdot_max_struct.max_xdot);
        xdot_max_struct.max_xdot_std = std(xdot_max_struct.max_xdot);
        xdot_max_struct.max_xdot_max = max(xdot_max_struct.max_xdot);
        xdot_max_struct.nTrials = nTrialsSubject;
        xdot_max_struct.SubjectNumber = subjectNumber;
        
        res = [res, xdot_max_struct];


end