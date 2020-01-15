function max_mean = maxMVC(mvc_dat)
% purpose: to find the maximum voluntary contraction of the subject using the sliding window method.
% This is done because taking the max-value from the EMG data would be very inaccurate. Instead, we will
% define a window that will pass through the data and largest mean of the window will be our max MVC 
% input: 
% mvc_dat = vector of maximum volume contraction data (single column of data) 
% output:
% max_mean = maximum of mvc data (a scalar)
% max 

data = length(mvc_dat); % vector of data 
n = 1500; % window size 
k = data - n; %creating the number of steps, the window will move through the dataone point(step) at a time
max_mean = 0;
% iterating through the steps of the window passing through the data 
for i = 1:k 
    
   % moving the window
   maximum = mean(mvc_dat(i:i+n-1,:));
   
   % reassigning max when necessary
   if maximum > max_mean
       max_mean(i) = maximum;
   end 
max_mean = max(max_mean);
end 