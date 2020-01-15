function res = GetZeroLevel(data)
    % GETZEROLEVEL gets an input of a single signal of 1xn samples, called
    % data.  It returns the zero level, which is defined as the mean of the
    % first samples.
    firstNPoints = 30;
    res = mean(data(1:firstNPoints));

end