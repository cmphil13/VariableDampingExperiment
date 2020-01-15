%% Plots for Phillips Fall19 FURI Poster 
%%
% Purpose of this code is to take in EMG data from trials variable of ImportData and represent
% the percent activation of the muscle compared to MVC in each of the damping 
% conditions: negative, positive, and variable. The data from each set of trials will be averaged
% and divided by the maximum voluntary contraction determined individually for each muscle
%% Initializing Variables and Clearing Residual Data 
clear emg_signals filt_emg_signals norm_emg maximum min_dat min_data filt_mvc_signals 
clear fore_pos bi_pos front_delt_pos rear_delt_pos long_tri_pos lat_tri_pos 
clear fore_neg bi_neg front_delt_neg rear_delt_neg long_tri_neg lat_tri_neg
clear fore_var bi_var front_delt_var rear_delt_var long_tri_var lat_tri_var

clear fore_avg_pos bi_avg_pos front_delt_avg_pos rear_delt_avg_pos long_tri_avg_pos lat_tri_avg_pos 
clear fore_avg_neg bi_avg_neg front_delt_avg_neg rear_delt_avg_neg long_tri_avg_neg lat_tri_avg_neg
clear fore_avg_var bi_avg_var front_delt_avg_var rear_delt_avg_var long_tri_avg_var lat_tri_avg_var

% loading MVC data 
mvc = load('MVC.txt');

x_axis_mvc = [1:length(mvc)]'; % axis for mvc 
mvc = mvc(:,1:6); %six because of the six different sensors, removing columns of zeros  

figure; plot(x_axis_mvc,mvc(:,1));
% initializing counts for later 
count_fore_var = 0;
count_bi_var = 0;
count_front_delt_var = 0;
count_rear_delt_var = 0;
count_long_tri_var = 0;
count_lat_tri_var = 0;

count_fore_neg = 0;
count_bi_neg = 0;
count_front_delt_neg = 0;
count_rear_delt_neg = 0;
count_long_tri_neg = 0;
count_lat_tri_neg = 0;

count_fore_pos = 0;
count_bi_pos = 0;
count_front_delt_pos = 0;
count_rear_delt_pos = 0;
count_long_tri_pos = 0;
count_lat_tri_pos = 0;

% Removing calibration trials 0,1,2,3 from both up and down
% trials_slim = trials;
% 
% for l = 1:length(trials_slim)
%     if trials_slim(l).GroupNumber == 0 || 1 || 2 || 3
%         trials_slim(l) = 0;
%     end
% end 


% Determining the minimum length of EMG data 
for k = 1:length(trials)
min_dat(k) = min(length(trials(k).Data.Emg));
end 

min_data = mink(min_dat,1); % this will help determine a uniform length to truncate the EMG data and also to plot against

%%
% iterating through number of trials 
for k = 190:380
   
% reading in emg data for the trial and trimming file to only include the data from the sensors columns [7:16] are all zeros 
emg_signals = trials(k).Data.Emg(1:min_data,1:6);

x_axis = 1:min_data; % to plot against data points 

for i = 1:6 % because there are six sensors 
    
% remove mean and take abs of mvc
mvc(:,i) = abs(mvc(:,i) - mean(mvc(:,i)));

% filter MVC data through low pass Butterworth filter 
filt_mvc_signals(:,i) = Butterworth_LPF(mvc(:,i), 1000, 1, 4); 
figure; plot(x_axis_mvc, filt_mvc_signals(:,1)); title('Filtered EMG data');xlabel('No. of Data Points');ylabel('Signal Amplitude');
figure; plot(x_axis_mvc, mvc(:,1)); title('Un-Filtered EMG data');xlabel('No. of Data Points');ylabel('Signal Amplitude');

% find max mvc 
% plot(x_axis_mvc, filt_mvc_signals(:,i)); % plot MVC
maximum =  max(filt_mvc_signals(:,i));

% remove mean and take abs of emg 
emg_signals(:,i) = abs(emg_signals(:,i) - mean(emg_signals(:,i)));

% filter EMG data through low pass Butterworth filter 
filt_emg_signals(:,i) = Butterworth_LPF(emg_signals(:,i), 1000, 20, 4);

% normalize EMG to MVC 
norm_emg(:,i) = filt_emg_signals(:,i)/maximum * 100; %percent activation
unfilt_norm(:,i) = emg_signals(:,i)/maximum * 100;

% Binning data based on type of damping, one if statement for each damping condition and one for each muscle 
if trials(k).DampingNumber == 3
    if i == 1
      count_fore_var = count_fore_var + 1;
      fore_var(:,count_fore_var) = norm_emg(:,i);
    elseif i == 2
      count_bi_var = count_bi_var + 1;
      bi_var(:,count_bi_var) = norm_emg(:,i);
    elseif i == 3 
      count_front_delt_var = count_front_delt_var + 1;
      front_delt_var(:,count_front_delt_var) = norm_emg(:,i);
    elseif i == 4 
      count_rear_delt_var = count_rear_delt_var + 1;
      rear_delt_var(:,count_rear_delt_var) = norm_emg(:,i);
    elseif i == 5
      count_long_tri_var = count_long_tri_var + 1; 
      long_tri_var(:,count_long_tri_var) = norm_emg(:,i);
    elseif i == 6
      count_lat_tri_var = count_lat_tri_var + 1;
      lat_tri_var(:,count_lat_tri_var) = norm_emg(:,i);
    end
end
if trials(k).DampingNumber == 2
    if i == 1
      count_fore_neg = count_fore_neg + 1;
      fore_neg(:,count_fore_neg) = norm_emg(:,i);
    elseif i == 2
      count_bi_neg = count_bi_neg + 1;
      bi_neg(:,count_bi_neg) = norm_emg(:,i);
    elseif i == 3 
      count_front_delt_neg = count_front_delt_neg + 1;
      front_delt_neg(:,count_front_delt_neg) = norm_emg(:,i);
    elseif i == 4 
      count_rear_delt_neg = count_rear_delt_neg + 1;
      rear_delt_neg(:,count_rear_delt_neg) = norm_emg(:,i);
    elseif i == 5
      count_long_tri_neg = count_long_tri_neg + 1; 
      long_tri_neg(:,count_long_tri_neg) = norm_emg(:,i);
    elseif i == 6
      count_lat_tri_neg = count_lat_tri_neg + 1;
      lat_tri_neg(:,count_lat_tri_neg) = norm_emg(:,i);
    end
end
if trials(k).DampingNumber == 1
    if i == 1
      count_fore_pos = count_fore_pos + 1;
      fore_pos(:,count_fore_pos) = norm_emg(:,i);
    elseif i == 2
      count_bi_pos = count_bi_pos + 1;
      bi_pos(:,count_bi_pos) = norm_emg(:,i);
    elseif i == 3 
      count_front_delt_pos = count_front_delt_pos + 1;
      front_delt_pos(:,count_front_delt_pos) = norm_emg(:,i);
    elseif i == 4 
      count_rear_delt_pos = count_rear_delt_pos + 1;
      rear_delt_pos(:,count_rear_delt_pos) = norm_emg(:,i);
    elseif i == 5
      count_long_tri_pos = count_long_tri_pos + 1; 
      long_tri_pos(:,count_long_tri_pos) = norm_emg(:,i);
    elseif i == 6
      count_lat_tri_pos = count_lat_tri_pos + 1;
      lat_tri_pos(:,count_lat_tri_pos) = norm_emg(:,i);
    end
end 
end
end


%% Getting Averages 
% % for j = 1:3 
% % averages.forearm.damping = [fore_avg_pos,fore_avg_neg,fore_avg_var]
% % averages.bicep.damping = [
% if mu
% for 
% averages.direction = ;
% averages.muscle = muscles(i);
% averages.damping =
% averages.damping.neg
% averages.damping.pos
% averages.damping.var
% averages.direction(1).muscle(2) = 'Bicep'
% averages.direction(1).muscle(3) = ''
% 
% averages.muscle.front_delt.damping.pos
% averages.muscle.rear_delt.damping
% averages.muscle.long_tri.damping
% averages.muscle.lat_tri.damping 
%     
    
fore_avg_var = mean(fore_var');
fore_avg_var = fore_avg_var';
fore_avg_pos = mean(fore_pos');
fore_avg_pos = fore_avg_pos';
fore_avg_neg = mean(fore_neg');
fore_avg_neg = fore_avg_neg';

bi_avg_var = mean(bi_var');
bi_avg_var = bi_avg_var';
bi_avg_pos = mean(bi_pos');
bi_avg_pos = bi_avg_pos';
bi_avg_neg = mean(bi_neg');
bi_avg_neg = bi_avg_neg';

front_delt_avg_var = mean(front_delt_var');
front_delt_avg_var = front_delt_avg_var';
front_delt_avg_pos = mean(front_delt_pos');
front_delt_avg_pos = front_delt_avg_pos';
front_delt_avg_neg = mean(front_delt_neg');
front_delt_avg_neg = front_delt_avg_neg';

rear_delt_avg_var = mean(rear_delt_var');
rear_delt_avg_var = rear_delt_avg_var';
rear_delt_avg_pos = mean(rear_delt_pos');
rear_delt_avg_pos = rear_delt_avg_pos';
rear_delt_avg_neg = mean(rear_delt_neg');
rear_delt_avg_neg = rear_delt_avg_neg';

long_tri_avg_var = mean(long_tri_var');
long_tri_avg_var = long_tri_avg_var';
long_tri_avg_pos = mean(long_tri_pos');
long_tri_avg_pos = long_tri_avg_pos';
long_tri_avg_neg = mean(long_tri_neg');
long_tri_avg_neg = long_tri_avg_neg';

lat_tri_avg_var = mean(lat_tri_var');
lat_tri_avg_var = lat_tri_avg_var';
lat_tri_avg_pos = mean(lat_tri_pos');
lat_tri_avg_pos = lat_tri_avg_pos';
lat_tri_avg_neg = mean(lat_tri_neg');
lat_tri_avg_neg = lat_tri_avg_neg';


% implement graph title based on muscle name
% for i =1:6 
% muscle_names = ["Forearm","Bicep","Front Delt","Rear Delt","Long. Tricep","Lat. Tricep"];
% graph_title = sprintf('Muscle Activation of %s',muscle_names(i)); 
figure;
plot(x_axis,fore_avg_pos,x_axis,fore_avg_neg,x_axis,fore_avg_var); title("Muscle Activation of Forearm for Different Damping Conditions"); xlabel("# of Data Points"); ylabel("Percent Activation");legend('Positive','Negative','Variable');
figure;
plot(x_axis,bi_avg_pos,x_axis,bi_avg_neg,x_axis,bi_avg_var); title("Muscle Activation of Bicep for Different Damping Conditions"); xlabel("# of Data Points"); ylabel("Percent Activation");legend('Positive','Negative','Variable');
figure;
plot(x_axis,front_delt_avg_pos,x_axis,front_delt_avg_neg,x_axis,front_delt_avg_var); title("Muscle Activation of Front Delt. for Different Damping Conditions"); xlabel("# of Data Points"); ylabel("Percent Activation");legend('Positive','Negative','Variable');
figure;
plot(x_axis,rear_delt_avg_pos,x_axis,rear_delt_avg_neg,x_axis,rear_delt_avg_var); title("Muscle Activation of Rear Delt. for Different Damping Conditions"); xlabel("# of Data Points"); ylabel("Percent Activation");legend('Positive','Negative','Variable');
figure;
plot(x_axis,long_tri_avg_pos,x_axis,long_tri_avg_neg,x_axis,long_tri_avg_var); title("Muscle Activation of Long. Tri. for Different Damping Conditions"); xlabel("# of Data Points"); ylabel("Percent Activation");legend('Positive','Negative','Variable');
figure;
plot(x_axis,lat_tri_avg_pos,x_axis,lat_tri_avg_neg,x_axis,lat_tri_avg_var); title("Muscle Activation of Lat. Tri. for Different Damping Conditions"); xlabel("# of Data Points"); ylabel("Percent Activation");legend('Positive','Negative','Variable');
