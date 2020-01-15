function res = CalcRepeatedAnova(data)
    % CALCREPEATEDANOVA takes in an array of data that nSubjects x
    % nConditions.  
    
    [nSubjects, nConditions] = size(data);
    
    % Calc Degrees of Freedom
    DFcond= nConditions - 1;
    DFsub = nSubjects - 1;
    DFerror = DFcond*DFsub;
    
    % Calc means
    subjectMeans = mean(data, 2);
    conditionMeans = mean(data, 1);
    grandMean = mean(subjectMeans);

    % Sum of squares error - condition
    SScond = nSubjects*(sum((conditionMeans - grandMean).^2));
    
    % Sum of squares error - within groups
    SSw = 0;
    for iCond = 1:nConditions
        SSw = SSw + sum((data(:,iCond) - conditionMeans(iCond)).^2);
    end
    
    % Sum of squares error - subjects
    SSsub = nConditions*(sum((subjectMeans - grandMean).^2));
    
    % Sum of squares error
    SSerror = SSw - SSsub;
    
    % Mean sum of squares - condition
    MScond = SScond/DFcond;
    
    % Mean sum of squares - error
    MSerror = SSerror/(DFsub*DFcond);
    
    % F statistic
    F = MScond/MSerror;
    
%     fprintf('Source\tSS\tdf\tMS\tF\n')
%     fprintf('Cond.\t%f\t%d\t%f\t%f\n', SScond, DFcond, MScond, F);
%     fprintf('Error\t%f\t%d\t%f\n\n', SSerror, DFerror, MSerror); 
    
    p = 1-fcdf(F, DFcond, DFerror);
    
%     fprintf('F(%d,%d) = %f\tp = %f\n', DFcond, DFerror, F, p);
    
    res = struct;
    res.F = F;
    res.DFcond = DFcond;
    res.DFerror = DFerror;
    res.p = p;
    
end

