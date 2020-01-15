function res = GetFivePercentRiseInd(data)
    
    data_temp = abs(data);
    zero_level = GetZeroLevel(data_temp);
    ss_level = GetSteadyStateLevel(data_temp);
    
    data_range = abs(ss_level-zero_level);
    

    fp = 0.05; % five percent
    inds = data_temp < (zero_level + fp*data_range);
    
    res = find(inds, 1, 'last');
    
end