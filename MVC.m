mvc = load('MVC_Tanner_28.txt');
% mvc = mvc(:,1:5); %six because of the six different sensors, removing columns of zeros  
% test_axis = length(mvc);

t_final_mvc = length(mvc(:,1))*1/2000; % getting time value for MVC
x_axis_mvc = [0.0005:0.0005:t_final_mvc]'; % creating vector for plot 
%% 
mvc = load('MVC.txt');
% mvc = mvc(:,1:5); %six because of the six different sensors, removing columns of zeros  
% test_axis = length(mvc);

t_final_mvc = length(mvc(:,1))*1/2000; % getting time value for MVC
x_axis_mvc = [0.0005:0.0005:t_final_mvc]'; % creating vector for plot 
%%
for i = 1:6
% remove mean and take abs of mvc
rect_mvc(:,i) = abs(mvc(:,i) - mean(mvc(:,i)));

% filter MVC data through low pass Butterworth filter 
filt_mvc_signals(:,i) = Butterworth_LPF(rect_mvc(:,i), 2000, 5, 4); 


% plots of mvc data
figure;
muscles = ["Forearm","Bicep","Front Delt","Rear Delt","Lateral Tricep","Longitudinal Tricep"];
title_mvc_1 = sprintf('MVC Signal Comparison for %s',muscles(i));
plot(x_axis_mvc, rect_mvc(:,i));xlabel('Time (s)');ylabel('Signal Amplitude (mV)'); 
hold on;
plot(x_axis_mvc, filt_mvc_signals(:,i),'r','linewidth',2); title(title_mvc_1);xlabel('Time (s)');ylabel('Signal Amplitude (mV)'); 

% legend('Filtered','Un-Filtered');
end

%% 
for i = 1:2
    
% remove mean and take abs of mvc
mvc(:,i) = abs(mvc(:,i) - mean(mvc(:,i)));

filt_mvc_signals(:,i) = Butterworth_LPF(mvc(:,i), 2000, 5, 4); 

muscles = ["Forearm","Bicep","Front Delt","Rear Delt","Lateral Tricep","Longitudinal Tricep"];
title_mvc_filt = sprintf('TFiltered MVC Signals for %s',muscles(i));
title_mvc = sprintf('TUnfiltered MVC Signals for %s',muscles(i));
figure;
plot(x_axis_mvc, filt_mvc_signals(:,i)); title(title_mvc_filt);xlabel('Time (s)');ylabel('Signal Amplitude (mV)');
% figure;
% plot(x_axis_mvc, mvc(:,1)); title(title_mvc);xlabel('Time (s)');ylabel('Signal Amplitude (mV)'); 
% legend('Filtered','Un-Filtered');
end 
