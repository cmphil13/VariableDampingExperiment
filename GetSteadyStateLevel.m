function res = GetSteadyStateLevel(data)
    % GETSTEADYSTATELEVEL gets an input of a single signal of 1xn samples, called
    % data.  It returns the steady state level, which is defined as the mean of the
    % last 100 samples.
    lastNPoints = 100;
    res = mean(data(end-lastNPoints:end));

end