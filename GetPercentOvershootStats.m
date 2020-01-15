function res = GetPercentOvershootStats(trials_wov, nSamples)
    % GETPERCENTOVERSHOOTSTATS takes in a the whole trial struct of trials
    % without violations and nSamples.  It returns a struct array of the
    % individual trial percent overshoot (po), max po, mean po, std po for
    % each subject.
    
    res = [];
    
    subjectNumbersAll = [trials_wov.SubjectNumber];
    subjectNumbers = unique(subjectNumbersAll);
    
    for i = 1:length(subjectNumbers)
        subjectNumber = subjectNumbers(i);
        subjectInds = subjectNumbersAll == subjectNumber;
        
        trials_subject = trials_wov(subjectInds);
        nTrialsSubject = length(trials_subject);
        x = zeros(nSamples, nTrialsSubject);
        for j = 1:nTrialsSubject
            x(:,j) = trials_subject(j).Data.EndEffPos_FromJA(1:nSamples);
        end
        
        trial_zero = mean(x(1:30,:), 1);
        trial_ss = mean(x(end-100:end,:), 1);
        trial_range = abs(trial_ss - trial_zero);
        trial_max = max(abs(x));
        
        po_struct = struct;       
        
        po_struct.po = (abs(trial_max - trial_zero)./trial_range - 1)*100;
        po_struct.po_mean = mean(po_struct.po);
        po_struct.po_std = std(po_struct.po);
        po_struct.po_max = max(po_struct.po);
        po_struct.nTrials = nTrialsSubject;
        po_struct.SubjectNumber = subjectNumber;
        
        res = [res, po_struct];
        
        
    end
    
    
    

    

    
end