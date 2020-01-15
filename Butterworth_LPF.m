% Butterworth code given by Dr. Lee, this is the framework for a low-pass filter  
function Signal_filtered = Butterworth_LPF(Signal, frequency, cut_frequency, order)
%input:
% Signal = data you want to filter, assumed as vector 
% frequency (kHz) = frequency of data you're inputting 
% cut_frequency = frequency you want to cut out of your data, closer to 0 is smoother, higher is more rough  
% order = order of the filter, higher order is more accurate and may be too small to visualize on plot, 4th order works well for EMG 
%output:
% Signal_filtered = filtered data  

% 5-7 Hz 

[b_low, a_low] = butter(order, ((cut_frequency) / (frequency/2)), 'low');
Signal_filtered = Signal;
for i = 1:size(Signal,2)
    Signal_filtered(:,i) = filtfilt(b_low, a_low, Signal(:,i));
end
end 


