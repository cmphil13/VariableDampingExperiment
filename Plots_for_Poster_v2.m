%% Plots for Phillips Fall19 FURI Poster 
%%
% Purpose of this code is to take in EMG data from trials variable of ImportData and represent
% the percent activation of the muscle compared to MVC in each of the damping 
% conditions: negative, positive, and variable. The data from each set of trials will be averaged
% and divided by the maximum voluntary contraction determined individually for each muscle

%% Order of Sections: 

% 1)  Import Subject.mat file and MVC
% 2)  Plot MVC 
% 3)  Intialize counts, clearing residual data, and defining min data range for EMG 
% 4)  Percent Activation Filtering and Parsing
% 5)  Movement Filtering and Parsing 
% 6)  Stability Filtering and Parsing 
% 7)  Mean Calculations 
% 8)  Percent Activation Plotting 
% 9)  Movement Bar Graphs for Average Muscle Activation
% 10) Stability Bar Graphs for Average Muscle Activation 

%% 1) Import Subject#Data.mat file and MVC
subjectMatFolder = './SubjectDataMats';

subjectMats = {
%                 'Subject20Data.mat', ...
%                 'Subject21Data.mat', ...
%                 'Subject22Data.mat', ...
%                 'Subject23Data.mat', ...
%                 'Subject24Data.mat', ...
%                 'Subject25Data.mat', ...
%                 'Subject26Data.mat', ...
%                 'Subject28Data.mat', ...
%                 'Subject29Data.mat', ...
%                 'Subject30Data.mat', ...
%                 'Subject31Data.mat', ...
                'Subject6Data.mat'
              };

clear trials
trials = [];
for i = 1:length(subjectMats)
    file = fullfile(subjectMatFolder, subjectMats{i});
    temp = load(file, 'trials');
    trials = [trials; temp.trials];
end

% Remember to change MVC file in 'VariableDampingAnalysis-master'
% loading MVC data 
mvc = load('MVC_6.txt');
mvc = mvc(:,1:6); %six because of the six different sensors, removing columns of zeros  
% test_axis = length(mvc);

t_final_mvc = length(mvc(:,1))*1/2000; % getting time value for MVC
x_axis_mvc = [0.0005:0.0005:t_final_mvc]'; % creating vector for plot 

%% 3) Initializing Variables and Clearing Residual Data 

% Clearing base signals 
clear emg_signals rect_mvc filt_emg_signals norm_emg maximum min_dat min_data filt_mvc_signals x_axis
clear fore_pos_4 bi_pos_4 front_delt_pos_4 rear_delt_pos_4 long_tri_pos_4 lat_tri_pos_4 
clear fore_neg_4 bi_neg_4 front_delt_neg_4 rear_delt_neg_4 long_tri_neg_4 lat_tri_neg_4
clear fore_var_4 bi_var_4 front_delt_var_4 rear_delt_var_4 long_tri_var_4 lat_tri_var_4

clear fore_pos_3 bi_pos_3 front_delt_pos_3 rear_delt_pos_3 long_tri_pos_3 lat_tri_pos_3 
clear fore_neg_3 bi_neg_3 front_delt_neg_3 rear_delt_neg_3 long_tri_neg_3 lat_tri_neg_3
clear fore_var_3 bi_var_3 front_delt_var_3 rear_delt_var_3 long_tri_var_3 lat_tri_var_3

clear fore_pos_2 bi_pos_2 front_delt_pos_2 rear_delt_pos_2 long_tri_pos_2 lat_tri_pos_2 
clear fore_neg_2 bi_neg_2 front_delt_neg_2 rear_delt_neg_2 long_tri_neg_2 lat_tri_neg_2
clear fore_var_2 bi_var_2 front_delt_var_2 rear_delt_var_2 long_tri_var_2 lat_tri_var_2

clear fore_pos_1 bi_pos_1 front_delt_pos_1 rear_delt_pos_1 long_tri_pos_1 lat_tri_pos_1 
clear fore_neg_1 bi_neg_1 front_delt_neg_1 rear_delt_neg_1 long_tri_neg_1 lat_tri_neg_1
clear fore_var_1 bi_var_1 front_delt_var_1 rear_delt_var_1 long_tri_var_1 lat_tri_var_1

% Clearing averages 
clear fore_avg_pos_4 bi_avg_pos_4 front_delt_avg_pos_4 rear_delt_avg_pos_4 long_tri_avg_pos_4 lat_tri_avg_pos_4 
clear fore_avg_neg_4 bi_avg_neg_4 front_delt_avg_neg_4 rear_delt_avg_neg_4 long_tri_avg_neg_4 lat_tri_avg_neg_4
clear fore_avg_var_4 bi_avg_var_4 front_delt_avg_var_4 rear_delt_avg_var_4 long_tri_avg_var_4 lat_tri_avg_var_4

clear fore_avg_pos_3 bi_avg_pos_3 front_delt_avg_pos_3 rear_delt_avg_pos_3 long_tri_avg_pos_3 lat_tri_avg_pos_3 
clear fore_avg_neg_3 bi_avg_neg_3 front_delt_avg_neg_3 rear_delt_avg_neg_3 long_tri_avg_neg_3 lat_tri_avg_neg_3
clear fore_avg_var_3 bi_avg_var_3 front_delt_avg_var_3 rear_delt_avg_var_3 long_tri_avg_var_3 lat_tri_avg_var_3

clear fore_avg_pos_2 bi_avg_pos_2 front_delt_avg_pos_2 rear_delt_avg_pos_2 long_tri_avg_pos_2 lat_tri_avg_pos_2 
clear fore_avg_neg_2 bi_avg_neg_2 front_delt_avg_neg_2 rear_delt_avg_neg_2 long_tri_avg_neg_2 lat_tri_avg_neg_2
clear fore_avg_var_2 bi_avg_var_2 front_delt_avg_var_2 rear_delt_avg_var_2 long_tri_avg_var_2 lat_tri_avg_var_2

clear fore_avg_pos_1 bi_avg_pos_1 front_delt_avg_pos_1 rear_delt_avg_pos_1 long_tri_avg_pos_1 lat_tri_avg_pos_1 
clear fore_avg_neg_1 bi_avg_neg_1 front_delt_avg_neg_1 rear_delt_avg_neg_1 long_tri_avg_neg_1 lat_tri_avg_neg_1
clear fore_avg_var_1 bi_avg_var_1 front_delt_avg_var_1 rear_delt_avg_var_1 long_tri_avg_var_1 lat_tri_avg_var_1

% initializing counts for binning data 
% upwards
count_fore_var_4 = 0;
count_bi_var_4 = 0;
count_front_delt_var_4 = 0;
count_rear_delt_var_4 = 0;
count_long_tri_var_4 = 0;
count_lat_tri_var_4 = 0;

count_fore_neg_4 = 0;
count_bi_neg_4 = 0;
count_front_delt_neg_4 = 0;
count_rear_delt_neg_4 = 0;
count_long_tri_neg_4 = 0;
count_lat_tri_neg_4 = 0;

count_fore_pos_4 = 0;
count_bi_pos_4 = 0;
count_front_delt_pos_4 = 0;
count_rear_delt_pos_4 = 0;
count_long_tri_pos_4 = 0;
count_lat_tri_pos_4 = 0;

%downwards
count_fore_var_3 = 0;
count_bi_var_3 = 0;
count_front_delt_var_3 = 0;
count_rear_delt_var_3 = 0;
count_long_tri_var_3 = 0;
count_lat_tri_var_3 = 0;

count_fore_neg_3 = 0;
count_bi_neg_3 = 0;
count_front_delt_neg_3 = 0;
count_rear_delt_neg_3 = 0;
count_long_tri_neg_3 = 0;
count_lat_tri_neg_3 = 0;

count_fore_pos_3 = 0;
count_bi_pos_3 = 0;
count_front_delt_pos_3 = 0;
count_rear_delt_pos_3 = 0;
count_long_tri_pos_3 = 0;
count_lat_tri_pos_3 = 0;

%right
count_fore_var_2 = 0;
count_bi_var_2 = 0;
count_front_delt_var_2 = 0;
count_rear_delt_var_2 = 0;
count_long_tri_var_2 = 0;
count_lat_tri_var_2 = 0;

count_fore_neg_2 = 0;
count_bi_neg_2 = 0;
count_front_delt_neg_2 = 0;
count_rear_delt_neg_2 = 0;
count_long_tri_neg_2 = 0;
count_lat_tri_neg_2 = 0;

count_fore_pos_2 = 0;
count_bi_pos_2 = 0;
count_front_delt_pos_2 = 0;
count_rear_delt_pos_2 = 0;
count_long_tri_pos_2 = 0;
count_lat_tri_pos_2 = 0;

%left 
count_fore_var_1 = 0;
count_bi_var_1 = 0;
count_front_delt_var_1 = 0;
count_rear_delt_var_1 = 0;
count_long_tri_var_1 = 0;
count_lat_tri_var_1 = 0;

count_fore_neg_1 = 0;
count_bi_neg_1 = 0;
count_front_delt_neg_1 = 0;
count_rear_delt_neg_1 = 0;
count_long_tri_neg_1 = 0;
count_lat_tri_neg_1 = 0;

count_fore_pos_1 = 0;
count_bi_pos_1 = 0;
count_front_delt_pos_1 = 0;
count_rear_delt_pos_1 = 0;
count_long_tri_pos_1 = 0;
count_lat_tri_pos_1 = 0;
% Whole subject 

% Determining the minimum length of EMG data 
for k = 1:380
min_dat(k) = min(length(trials(k).Data.Emg));
end 
 min_data = mink(min_dat,1);

% For Subject 6 ML 
min_data =(min_dat(191:380));
min_data_p =mink(min_data,1);% this will help determine a uniform length to truncate the EMG data and also to plot against

%% 4)  Percent Activation Filtering Data and Seperating to Vectors
% iterating through number of trials 
for k = 1:190
   
% reading in emg data for the trial and trimming file to only include the data from the sensors columns [7:16] are all zeros 
emg_signals = trials(k).Data.Emg(1:min_data,1:6);
t_final = min_data*1/2000; 
x_axis = [0.0005:0.0005:t_final]'; % since sampling frequency is 2000 Hz the spacing between points is 1/2000 = 0.0005 s 

if trials(k).GroupNumber ~=2 
if trials(k).GroupNumber ~=3

for i = 1:6 % because there are six sensors 
    
% remove mean and take abs of mvc
rect_mvc(:,i) = abs(mvc(:,i) - mean(mvc(:,i)));

% filter MVC data through low pass Butterworth filter 
filt_mvc_signals(:,i) = Butterworth_LPF(rect_mvc(:,i), 2000, 5, 4); 

% Finding MVC for each muscle 
max_val =  maxk(filt_mvc_signals(:,i),10000);
maximum(i) = max_val(10000);

% remove mean and take abs of emg 
emg_signals(:,i) = abs(emg_signals(:,i) - mean(emg_signals(:,i)));

% filter EMG data through low pass Butterworth filter 
filt_emg_signals(:,i) = Butterworth_LPF(emg_signals(:,i), 2000, 5, 4);

% normalize EMG to MVC 
norm_emg(:,i) = filt_emg_signals(:,i)/maximum(i) * 100; %percent activation
% Binning data based on type of damping, one if statement for each damping condition and one for each muscle 
if trials(k).TargetDirNum == 4

if trials(k).DampingNumber == 3
    if i == 1
      count_fore_var_4 = count_fore_var_4 + 1;
      fore_var_4(:,count_fore_var_4) = norm_emg(:,i);
    elseif i == 2
      count_bi_var_4 = count_bi_var_4 + 1;
      bi_var_4(:,count_bi_var_4) = norm_emg(:,i);
    elseif i == 3 
      count_front_delt_var_4 = count_front_delt_var_4 + 1;
      front_delt_var_4(:,count_front_delt_var_4) = norm_emg(:,i);
    elseif i == 4 
      count_rear_delt_var_4 = count_rear_delt_var_4 + 1;
      rear_delt_var_4(:,count_rear_delt_var_4) = norm_emg(:,i);
    elseif i == 5
      count_long_tri_var_4 = count_long_tri_var_4 + 1; 
      long_tri_var_4(:,count_long_tri_var_4) = norm_emg(:,i);
    elseif i == 6
      count_lat_tri_var_4 = count_lat_tri_var_4 + 1;
      lat_tri_var_4(:,count_lat_tri_var_4) = norm_emg(:,i);
    end
end

if trials(k).DampingNumber == 2
    if i == 1
      count_fore_neg_4 = count_fore_neg_4 + 1;
      fore_neg_4(:,count_fore_neg_4) = norm_emg(:,i);
    elseif i == 2
      count_bi_neg_4 = count_bi_neg_4 + 1;
      bi_neg_4(:,count_bi_neg_4) = norm_emg(:,i);
    elseif i == 3 
      count_front_delt_neg_4 = count_front_delt_neg_4 + 1;
      front_delt_neg_4(:,count_front_delt_neg_4) = norm_emg(:,i);
    elseif i == 4 
      count_rear_delt_neg_4 = count_rear_delt_neg_4 + 1;
      rear_delt_neg_4(:,count_rear_delt_neg_4) = norm_emg(:,i);
    elseif i == 5
      count_long_tri_neg_4 = count_long_tri_neg_4 + 1; 
      long_tri_neg_4(:,count_long_tri_neg_4) = norm_emg(:,i);
    elseif i == 6
      count_lat_tri_neg_4 = count_lat_tri_neg_4 + 1;
      lat_tri_neg_4(:,count_lat_tri_neg_4) = norm_emg(:,i);
    end
end
if trials(k).DampingNumber == 1
    if i == 1
      count_fore_pos_4 = count_fore_pos_4 + 1;
      fore_pos_4(:,count_fore_pos_4) = norm_emg(:,i);
    elseif i == 2
      count_bi_pos_4 = count_bi_pos_4 + 1;
      bi_pos_4(:,count_bi_pos_4) = norm_emg(:,i);
    elseif i == 3 
      count_front_delt_pos_4 = count_front_delt_pos_4 + 1;
      front_delt_pos_4(:,count_front_delt_pos_4) = norm_emg(:,i);
    elseif i == 4 
      count_rear_delt_pos_4 = count_rear_delt_pos_4 + 1;
      rear_delt_pos_4(:,count_rear_delt_pos_4) = norm_emg(:,i);
    elseif i == 5
      count_long_tri_pos_4 = count_long_tri_pos_4 + 1; 
      long_tri_pos_4(:,count_long_tri_pos_4) = norm_emg(:,i);
    elseif i == 6
      count_lat_tri_pos_4 = count_lat_tri_pos_4 + 1;
      lat_tri_pos_4(:,count_lat_tri_pos_4) = norm_emg(:,i);
    end
end 
end

if trials(k).TargetDirNum == 3

if trials(k).DampingNumber == 3
    if i == 1
      count_fore_var_3 = count_fore_var_3 + 1;
      fore_var_3(:,count_fore_var_3) = norm_emg(:,i);
    elseif i == 2
      count_bi_var_3 = count_bi_var_3 + 1;
      bi_var_3(:,count_bi_var_3) = norm_emg(:,i);
    elseif i == 3 
      count_front_delt_var_3 = count_front_delt_var_3 + 1;
      front_delt_var_3(:,count_front_delt_var_3) = norm_emg(:,i);
    elseif i == 4 
      count_rear_delt_var_3 = count_rear_delt_var_3 + 1;
      rear_delt_var_3(:,count_rear_delt_var_3) = norm_emg(:,i);
    elseif i == 5
      count_long_tri_var_3 = count_long_tri_var_3 + 1; 
      long_tri_var_3(:,count_long_tri_var_3) = norm_emg(:,i);
    elseif i == 6
      count_lat_tri_var_3 = count_lat_tri_var_3 + 1;
      lat_tri_var_3(:,count_lat_tri_var_3) = norm_emg(:,i);
    end
end

if trials(k).DampingNumber == 2
    if i == 1
      count_fore_neg_3 = count_fore_neg_3 + 1;
      fore_neg_3(:,count_fore_neg_3) = norm_emg(:,i);
    elseif i == 2
      count_bi_neg_3 = count_bi_neg_3 + 1;
      bi_neg_3(:,count_bi_neg_3) = norm_emg(:,i);
    elseif i == 3 
      count_front_delt_neg_3 = count_front_delt_neg_3 + 1;
      front_delt_neg_3(:,count_front_delt_neg_3) = norm_emg(:,i);
    elseif i == 4 
      count_rear_delt_neg_3 = count_rear_delt_neg_3 + 1;
      rear_delt_neg_3(:,count_rear_delt_neg_3) = norm_emg(:,i);
    elseif i == 5
      count_long_tri_neg_3 = count_long_tri_neg_3 + 1; 
      long_tri_neg_3(:,count_long_tri_neg_3) = norm_emg(:,i);
    elseif i == 6
      count_lat_tri_neg_3 = count_lat_tri_neg_3 + 1;
      lat_tri_neg_3(:,count_lat_tri_neg_3) = norm_emg(:,i);
    end
end
if trials(k).DampingNumber == 1
    if i == 1
      count_fore_pos_3 = count_fore_pos_3 + 1;
      fore_pos_3(:,count_fore_pos_3) = norm_emg(:,i);
    elseif i == 2
      count_bi_pos_3 =  count_bi_pos_3 + 1;
      bi_pos_3(:,count_bi_pos_3) = norm_emg(:,i);
    elseif i == 3 
      count_front_delt_pos_3 = count_front_delt_pos_3 + 1;
      front_delt_pos_3(:,count_front_delt_pos_3) = norm_emg(:,i);
    elseif i == 4 
      count_rear_delt_pos_3 = count_rear_delt_pos_3 + 1;
      rear_delt_pos_3(:,count_rear_delt_pos_3) = norm_emg(:,i);
    elseif i == 5
      count_long_tri_pos_3 = count_long_tri_pos_3 + 1; 
      long_tri_pos_3(:,count_long_tri_pos_3) = norm_emg(:,i);
    elseif i == 6
      count_lat_tri_pos_3 = count_lat_tri_pos_3 + 1;
      lat_tri_pos_3(:,count_lat_tri_pos_3) = norm_emg(:,i);
    end
end 
end

if trials(k).TargetDirNum == 2
if trials(k).DampingNumber == 3
    if i == 1
      count_fore_var_2 = count_fore_var_2 + 1;
      fore_var_2(:,count_fore_var_2) = norm_emg(:,i);
    elseif i == 2
      count_bi_var_2 = count_bi_var_2 + 1;
      bi_var_2(:,count_bi_var_2) = norm_emg(:,i);
    elseif i == 3 
      count_front_delt_var_2 = count_front_delt_var_2 + 1;
      front_delt_var_2(:,count_front_delt_var_2) = norm_emg(:,i);
    elseif i == 4 
      count_rear_delt_var_2 = count_rear_delt_var_2 + 1;
      rear_delt_var_2(:,count_rear_delt_var_2) = norm_emg(:,i);
    elseif i == 5
      count_long_tri_var_2 = count_long_tri_var_2 + 1; 
      long_tri_var_2(:,count_long_tri_var_2) = norm_emg(:,i);
    elseif i == 6
      count_lat_tri_var_2 = count_lat_tri_var_2 + 1;
      lat_tri_var_2(:,count_lat_tri_var_2) = norm_emg(:,i);
    end
end

if trials(k).DampingNumber == 2
    if i == 1
      count_fore_neg_2 = count_fore_neg_2 + 1;
      fore_neg_2(:,count_fore_neg_2) = norm_emg(:,i);
    elseif i == 2
      count_bi_neg_2 = count_bi_neg_2 + 1;
      bi_neg_2(:,count_bi_neg_2) = norm_emg(:,i);
    elseif i == 3 
      count_front_delt_neg_2 = count_front_delt_neg_2 + 1;
      front_delt_neg_2(:,count_front_delt_neg_2) = norm_emg(:,i);
    elseif i == 4 
      count_rear_delt_neg_2 = count_rear_delt_neg_2 + 1;
      rear_delt_neg_2(:,count_rear_delt_neg_2) = norm_emg(:,i);
    elseif i == 5
      count_long_tri_neg_2 = count_long_tri_neg_2 + 1; 
      long_tri_neg_2(:,count_long_tri_neg_2) = norm_emg(:,i);
    elseif i == 6
      count_lat_tri_neg_2 = count_lat_tri_neg_2 + 1;
      lat_tri_neg_2(:,count_lat_tri_neg_2) = norm_emg(:,i);
    end
end
if trials(k).DampingNumber == 1
    if i == 1
      count_fore_pos_2 = count_fore_pos_2 + 1;
      fore_pos_2(:,count_fore_pos_2) = norm_emg(:,i);
    elseif i == 2
      count_bi_pos_2 =  count_bi_pos_2 + 1;
      bi_pos_2(:,count_bi_pos_2) = norm_emg(:,i);
    elseif i == 3 
      count_front_delt_pos_2 = count_front_delt_pos_2 + 1;
      front_delt_pos_2(:,count_front_delt_pos_2) = norm_emg(:,i);
    elseif i == 4 
      count_rear_delt_pos_2 = count_rear_delt_pos_2 + 1;
      rear_delt_pos_2(:,count_rear_delt_pos_2) = norm_emg(:,i);
    elseif i == 5
      count_long_tri_pos_2 = count_long_tri_pos_2 + 1; 
      long_tri_pos_2(:,count_long_tri_pos_2) = norm_emg(:,i);
    elseif i == 6
      count_lat_tri_pos_2 = count_lat_tri_pos_2 + 1;
      lat_tri_pos_2(:,count_lat_tri_pos_2) = norm_emg(:,i);
    end
end 
end

if trials(k).TargetDirNum == 1
if trials(k).DampingNumber == 3
    if i == 1
      count_fore_var_1 = count_fore_var_1 + 1;
      fore_var_1(:,count_fore_var_1) = norm_emg(:,i);
    elseif i == 2
      count_bi_var_1 = count_bi_var_1 + 1;
      bi_var_1(:,count_bi_var_1) = norm_emg(:,i);
    elseif i == 3 
      count_front_delt_var_1 = count_front_delt_var_1 + 1;
      front_delt_var_1(:,count_front_delt_var_1) = norm_emg(:,i);
    elseif i == 4 
      count_rear_delt_var_1 = count_rear_delt_var_1 + 1;
      rear_delt_var_1(:,count_rear_delt_var_1) = norm_emg(:,i);
    elseif i == 5
      count_long_tri_var_1 = count_long_tri_var_1 + 1; 
      long_tri_var_1(:,count_long_tri_var_1) = norm_emg(:,i);
    elseif i == 6
      count_lat_tri_var_1 = count_lat_tri_var_1 + 1;
      lat_tri_var_1(:,count_lat_tri_var_1) = norm_emg(:,i);
    end
end

if trials(k).DampingNumber == 2
    if i == 1
      count_fore_neg_1 = count_fore_neg_1 + 1;
      fore_neg_1(:,count_fore_neg_1) = norm_emg(:,i);
    elseif i == 2
      count_bi_neg_1 = count_bi_neg_1 + 1;
      bi_neg_1(:,count_bi_neg_1) = norm_emg(:,i);
    elseif i == 3 
      count_front_delt_neg_1 = count_front_delt_neg_1 + 1;
      front_delt_neg_1(:,count_front_delt_neg_1) = norm_emg(:,i);
    elseif i == 4 
      count_rear_delt_neg_1 = count_rear_delt_neg_1 + 1;
      rear_delt_neg_1(:,count_rear_delt_neg_1) = norm_emg(:,i);
    elseif i == 5
      count_long_tri_neg_1 = count_long_tri_neg_1 + 1; 
      long_tri_neg_1(:,count_long_tri_neg_1) = norm_emg(:,i);
    elseif i == 6
      count_lat_tri_neg_1 = count_lat_tri_neg_1 + 1;
      lat_tri_neg_1(:,count_lat_tri_neg_1) = norm_emg(:,i);
    end
end
if trials(k).DampingNumber == 1
    if i == 1
      count_fore_pos_1 = count_fore_pos_1 + 1;
      fore_pos_1(:,count_fore_pos_1) = norm_emg(:,i);
    elseif i == 2
      count_bi_pos_1 =  count_bi_pos_1 + 1;
      bi_pos_1(:,count_bi_pos_1) = norm_emg(:,i);
    elseif i == 3 
      count_front_delt_pos_1 = count_front_delt_pos_1 + 1;
      front_delt_pos_1(:,count_front_delt_pos_1) = norm_emg(:,i);
    elseif i == 4 
      count_rear_delt_pos_1 = count_rear_delt_pos_1 + 1;
      rear_delt_pos_1(:,count_rear_delt_pos_1) = norm_emg(:,i);
    elseif i == 5
      count_long_tri_pos_1 = count_long_tri_pos_1 + 1; 
      long_tri_pos_1(:,count_long_tri_pos_1) = norm_emg(:,i);
    elseif i == 6
      count_lat_tri_pos_1 = count_lat_tri_pos_1 + 1;
      lat_tri_pos_1(:,count_lat_tri_pos_1) = norm_emg(:,i);
    end
end 
end

end 
end
end
end

%% 5) Filtering and Parsing for Movement Section 
% Filtering Data and Seperating to Vectors 
% iterating through number of trials 
for k = 191:380
   
% reading in emg data for the trial and trimming file to only include the data from the sensors columns [7:16] are all zeros 
emg_signals = trials(k).Data.Emg(500:1400,1:6); %truncating data to 0.25 and 0.75 
x_axis = [0.25:0.0005:0.70]'; % since sampling frequency is 2000 Hz the spacing between points is 1/2000 = 0.0005 s 

 
if trials(k).GroupNumber ~=2 
if trials(k).GroupNumber ~=3

for i = 1:6 % because there are six sensors 
    
% remove mean and take abs of mvc
rect_mvc(:,i) = abs(mvc(:,i) - mean(mvc(:,i)));

% filter MVC data through low pass Butterworth filter 
filt_mvc_signals(:,i) = Butterworth_LPF(rect_mvc(:,i), 2000, 5, 4); 

% Finding MVC for each muscle 
max_val =  maxk(filt_mvc_signals(:,i),20000);
maximum(i) = max_val(20000);

% remove mean and take abs of emg 
emg_signals(:,i) = abs(emg_signals(:,i) - mean(emg_signals(:,i)));

% filter EMG data through low pass Butterworth filter 
filt_emg_signals(:,i) = Butterworth_LPF(emg_signals(:,i), 2000, 5, 4);

% normalize EMG to MVC 
norm_emg(:,i) = filt_emg_signals(:,i)/maximum(i) * 100; %percent activation
% Binning data based on type of damping, one if statement for each damping condition and one for each muscle 
if trials(k).TargetDirNum == 4

if trials(k).DampingNumber == 3
    if i == 1
      count_fore_var_4 = count_fore_var_4 + 1;
      fore_var_4(:,count_fore_var_4) = norm_emg(:,i);
    elseif i == 2
      count_bi_var_4 = count_bi_var_4 + 1;
      bi_var_4(:,count_bi_var_4) = norm_emg(:,i);
    elseif i == 3 
      count_front_delt_var_4 = count_front_delt_var_4 + 1;
      front_delt_var_4(:,count_front_delt_var_4) = norm_emg(:,i);
    elseif i == 4 
      count_rear_delt_var_4 = count_rear_delt_var_4 + 1;
      rear_delt_var_4(:,count_rear_delt_var_4) = norm_emg(:,i);
    elseif i == 5
      count_long_tri_var_4 = count_long_tri_var_4 + 1; 
      long_tri_var_4(:,count_long_tri_var_4) = norm_emg(:,i);
    elseif i == 6
      count_lat_tri_var_4 = count_lat_tri_var_4 + 1;
      lat_tri_var_4(:,count_lat_tri_var_4) = norm_emg(:,i);
    end
end

if trials(k).DampingNumber == 2
    if i == 1
      count_fore_neg_4 = count_fore_neg_4 + 1;
      fore_neg_4(:,count_fore_neg_4) = norm_emg(:,i);
    elseif i == 2
      count_bi_neg_4 = count_bi_neg_4 + 1;
      bi_neg_4(:,count_bi_neg_4) = norm_emg(:,i);
    elseif i == 3 
      count_front_delt_neg_4 = count_front_delt_neg_4 + 1;
      front_delt_neg_4(:,count_front_delt_neg_4) = norm_emg(:,i);
    elseif i == 4 
      count_rear_delt_neg_4 = count_rear_delt_neg_4 + 1;
      rear_delt_neg_4(:,count_rear_delt_neg_4) = norm_emg(:,i);
    elseif i == 5
      count_long_tri_neg_4 = count_long_tri_neg_4 + 1; 
      long_tri_neg_4(:,count_long_tri_neg_4) = norm_emg(:,i);
    elseif i == 6
      count_lat_tri_neg_4 = count_lat_tri_neg_4 + 1;
      lat_tri_neg_4(:,count_lat_tri_neg_4) = norm_emg(:,i);
    end
end
if trials(k).DampingNumber == 1
    if i == 1
      count_fore_pos_4 = count_fore_pos_4 + 1;
      fore_pos_4(:,count_fore_pos_4) = norm_emg(:,i);
    elseif i == 2
      count_bi_pos_4 = count_bi_pos_4 + 1;
      bi_pos_4(:,count_bi_pos_4) = norm_emg(:,i);
    elseif i == 3 
      count_front_delt_pos_4 = count_front_delt_pos_4 + 1;
      front_delt_pos_4(:,count_front_delt_pos_4) = norm_emg(:,i);
    elseif i == 4 
      count_rear_delt_pos_4 = count_rear_delt_pos_4 + 1;
      rear_delt_pos_4(:,count_rear_delt_pos_4) = norm_emg(:,i);
    elseif i == 5
      count_long_tri_pos_4 = count_long_tri_pos_4 + 1; 
      long_tri_pos_4(:,count_long_tri_pos_4) = norm_emg(:,i);
    elseif i == 6
      count_lat_tri_pos_4 = count_lat_tri_pos_4 + 1;
      lat_tri_pos_4(:,count_lat_tri_pos_4) = norm_emg(:,i);
    end
end 
end

if trials(k).TargetDirNum == 3

if trials(k).DampingNumber == 3
    if i == 1
      count_fore_var_3 = count_fore_var_3 + 1;
      fore_var_3(:,count_fore_var_3) = norm_emg(:,i);
    elseif i == 2
      count_bi_var_3 = count_bi_var_3 + 1;
      bi_var_3(:,count_bi_var_3) = norm_emg(:,i);
    elseif i == 3 
      count_front_delt_var_3 = count_front_delt_var_3 + 1;
      front_delt_var_3(:,count_front_delt_var_3) = norm_emg(:,i);
    elseif i == 4 
      count_rear_delt_var_3 = count_rear_delt_var_3 + 1;
      rear_delt_var_3(:,count_rear_delt_var_3) = norm_emg(:,i);
    elseif i == 5
      count_long_tri_var_3 = count_long_tri_var_3 + 1; 
      long_tri_var_3(:,count_long_tri_var_3) = norm_emg(:,i);
    elseif i == 6
      count_lat_tri_var_3 = count_lat_tri_var_3 + 1;
      lat_tri_var_3(:,count_lat_tri_var_3) = norm_emg(:,i);
    end
end

if trials(k).DampingNumber == 2
    if i == 1
      count_fore_neg_3 = count_fore_neg_3 + 1;
      fore_neg_3(:,count_fore_neg_3) = norm_emg(:,i);
    elseif i == 2
      count_bi_neg_3 = count_bi_neg_3 + 1;
      bi_neg_3(:,count_bi_neg_3) = norm_emg(:,i);
    elseif i == 3 
      count_front_delt_neg_3 = count_front_delt_neg_3 + 1;
      front_delt_neg_3(:,count_front_delt_neg_3) = norm_emg(:,i);
    elseif i == 4 
      count_rear_delt_neg_3 = count_rear_delt_neg_3 + 1;
      rear_delt_neg_3(:,count_rear_delt_neg_3) = norm_emg(:,i);
    elseif i == 5
      count_long_tri_neg_3 = count_long_tri_neg_3 + 1; 
      long_tri_neg_3(:,count_long_tri_neg_3) = norm_emg(:,i);
    elseif i == 6
      count_lat_tri_neg_3 = count_lat_tri_neg_3 + 1;
      lat_tri_neg_3(:,count_lat_tri_neg_3) = norm_emg(:,i);
    end
end
if trials(k).DampingNumber == 1
    if i == 1
      count_fore_pos_3 = count_fore_pos_3 + 1;
      fore_pos_3(:,count_fore_pos_3) = norm_emg(:,i);
    elseif i == 2
      count_bi_pos_3 =  count_bi_pos_3 + 1;
      bi_pos_3(:,count_bi_pos_3) = norm_emg(:,i);
    elseif i == 3 
      count_front_delt_pos_3 = count_front_delt_pos_3 + 1;
      front_delt_pos_3(:,count_front_delt_pos_3) = norm_emg(:,i);
    elseif i == 4 
      count_rear_delt_pos_3 = count_rear_delt_pos_3 + 1;
      rear_delt_pos_3(:,count_rear_delt_pos_3) = norm_emg(:,i);
    elseif i == 5
      count_long_tri_pos_3 = count_long_tri_pos_3 + 1; 
      long_tri_pos_3(:,count_long_tri_pos_3) = norm_emg(:,i);
    elseif i == 6
      count_lat_tri_pos_3 = count_lat_tri_pos_3 + 1;
      lat_tri_pos_3(:,count_lat_tri_pos_3) = norm_emg(:,i);
    end
end 
end

if trials(k).TargetDirNum == 2
if trials(k).DampingNumber == 3
    if i == 1
      count_fore_var_2 = count_fore_var_2 + 1;
      fore_var_2(:,count_fore_var_2) = norm_emg(:,i);
    elseif i == 2
      count_bi_var_2 = count_bi_var_2 + 1;
      bi_var_2(:,count_bi_var_2) = norm_emg(:,i);
    elseif i == 3 
      count_front_delt_var_2 = count_front_delt_var_2 + 1;
      front_delt_var_2(:,count_front_delt_var_2) = norm_emg(:,i);
    elseif i == 4 
      count_rear_delt_var_2 = count_rear_delt_var_2 + 1;
      rear_delt_var_2(:,count_rear_delt_var_2) = norm_emg(:,i);
    elseif i == 5
      count_long_tri_var_2 = count_long_tri_var_2 + 1; 
      long_tri_var_2(:,count_long_tri_var_2) = norm_emg(:,i);
    elseif i == 6
      count_lat_tri_var_2 = count_lat_tri_var_2 + 1;
      lat_tri_var_2(:,count_lat_tri_var_2) = norm_emg(:,i);
    end
end

if trials(k).DampingNumber == 2
    if i == 1
      count_fore_neg_2 = count_fore_neg_2 + 1;
      fore_neg_2(:,count_fore_neg_2) = norm_emg(:,i);
    elseif i == 2
      count_bi_neg_2 = count_bi_neg_2 + 1;
      bi_neg_2(:,count_bi_neg_2) = norm_emg(:,i);
    elseif i == 3 
      count_front_delt_neg_2 = count_front_delt_neg_2 + 1;
      front_delt_neg_2(:,count_front_delt_neg_2) = norm_emg(:,i);
    elseif i == 4 
      count_rear_delt_neg_2 = count_rear_delt_neg_2 + 1;
      rear_delt_neg_2(:,count_rear_delt_neg_2) = norm_emg(:,i);
    elseif i == 5
      count_long_tri_neg_2 = count_long_tri_neg_2 + 1; 
      long_tri_neg_2(:,count_long_tri_neg_2) = norm_emg(:,i);
    elseif i == 6
      count_lat_tri_neg_2 = count_lat_tri_neg_2 + 1;
      lat_tri_neg_2(:,count_lat_tri_neg_2) = norm_emg(:,i);
    end
end
if trials(k).DampingNumber == 1
    if i == 1
      count_fore_pos_2 = count_fore_pos_2 + 1;
      fore_pos_2(:,count_fore_pos_2) = norm_emg(:,i);
    elseif i == 2
      count_bi_pos_2 =  count_bi_pos_2 + 1;
      bi_pos_2(:,count_bi_pos_2) = norm_emg(:,i);
    elseif i == 3 
      count_front_delt_pos_2 = count_front_delt_pos_2 + 1;
      front_delt_pos_2(:,count_front_delt_pos_2) = norm_emg(:,i);
    elseif i == 4 
      count_rear_delt_pos_2 = count_rear_delt_pos_2 + 1;
      rear_delt_pos_2(:,count_rear_delt_pos_2) = norm_emg(:,i);
    elseif i == 5
      count_long_tri_pos_2 = count_long_tri_pos_2 + 1; 
      long_tri_pos_2(:,count_long_tri_pos_2) = norm_emg(:,i);
    elseif i == 6
      count_lat_tri_pos_2 = count_lat_tri_pos_2 + 1;
      lat_tri_pos_2(:,count_lat_tri_pos_2) = norm_emg(:,i);
    end
end 
end

if trials(k).TargetDirNum == 1
if trials(k).DampingNumber == 3
    if i == 1
      count_fore_var_1 = count_fore_var_1 + 1;
      fore_var_1(:,count_fore_var_1) = norm_emg(:,i);
    elseif i == 2
      count_bi_var_1 = count_bi_var_1 + 1;
      bi_var_1(:,count_bi_var_1) = norm_emg(:,i);
    elseif i == 3 
      count_front_delt_var_1 = count_front_delt_var_1 + 1;
      front_delt_var_1(:,count_front_delt_var_1) = norm_emg(:,i);
    elseif i == 4 
      count_rear_delt_var_1 = count_rear_delt_var_1 + 1;
      rear_delt_var_1(:,count_rear_delt_var_1) = norm_emg(:,i);
    elseif i == 5
      count_long_tri_var_1 = count_long_tri_var_1 + 1; 
      long_tri_var_1(:,count_long_tri_var_1) = norm_emg(:,i);
    elseif i == 6
      count_lat_tri_var_1 = count_lat_tri_var_1 + 1;
      lat_tri_var_1(:,count_lat_tri_var_1) = norm_emg(:,i);
    end
end

if trials(k).DampingNumber == 2
    if i == 1
      count_fore_neg_1 = count_fore_neg_1 + 1;
      fore_neg_1(:,count_fore_neg_1) = norm_emg(:,i);
    elseif i == 2
      count_bi_neg_1 = count_bi_neg_1 + 1;
      bi_neg_1(:,count_bi_neg_1) = norm_emg(:,i);
    elseif i == 3 
      count_front_delt_neg_1 = count_front_delt_neg_1 + 1;
      front_delt_neg_1(:,count_front_delt_neg_1) = norm_emg(:,i);
    elseif i == 4 
      count_rear_delt_neg_1 = count_rear_delt_neg_1 + 1;
      rear_delt_neg_1(:,count_rear_delt_neg_1) = norm_emg(:,i);
    elseif i == 5
      count_long_tri_neg_1 = count_long_tri_neg_1 + 1; 
      long_tri_neg_1(:,count_long_tri_neg_1) = norm_emg(:,i);
    elseif i == 6
      count_lat_tri_neg_1 = count_lat_tri_neg_1 + 1;
      lat_tri_neg_1(:,count_lat_tri_neg_1) = norm_emg(:,i);
    end
end
if trials(k).DampingNumber == 1
    if i == 1
      count_fore_pos_1 = count_fore_pos_1 + 1;
      fore_pos_1(:,count_fore_pos_1) = norm_emg(:,i);
    elseif i == 2
      count_bi_pos_1 =  count_bi_pos_1 + 1;
      bi_pos_1(:,count_bi_pos_1) = norm_emg(:,i);
    elseif i == 3 
      count_front_delt_pos_1 = count_front_delt_pos_1 + 1;
      front_delt_pos_1(:,count_front_delt_pos_1) = norm_emg(:,i);
    elseif i == 4 
      count_rear_delt_pos_1 = count_rear_delt_pos_1 + 1;
      rear_delt_pos_1(:,count_rear_delt_pos_1) = norm_emg(:,i);
    elseif i == 5
      count_long_tri_pos_1 = count_long_tri_pos_1 + 1; 
      long_tri_pos_1(:,count_long_tri_pos_1) = norm_emg(:,i);
    elseif i == 6
      count_lat_tri_pos_1 = count_lat_tri_pos_1 + 1;
      lat_tri_pos_1(:,count_lat_tri_pos_1) = norm_emg(:,i);
    end
end 
end

end 
end
end
end

%% 6) Filtering and Parsing for Stability Section 
% Filtering Data and Seperating to Vectors 
% iterating through number of trials 
for k = 1:190
   
% reading in emg data for the trial and trimming file to only include the data from the sensors columns [7:16] are all zeros 
emg_signals = trials(k).Data.Emg(2000:2500,1:6); %truncating data to 0.25 and 0.75 
x_axis = [1:0.0005:1.5]'; % since sampling frequency is 2000 Hz the spacing between points is 1/2000 = 0.0005 s 

 
if trials(k).GroupNumber ~=2 
if trials(k).GroupNumber ~=3

for i = 1:6 % because there are six sensors 
    
% remove mean and take abs of mvc
rect_mvc(:,i) = abs(mvc(:,i) - mean(mvc(:,i)));

% filter MVC data through low pass Butterworth filter 
filt_mvc_signals(:,i) = Butterworth_LPF(rect_mvc(:,i), 2000, 5, 4); 

% Finding MVC for each muscle 
max_val =  maxk(filt_mvc_signals(:,i),30000);
maximum(i) = max_val(30000);

% remove mean and take abs of emg 
emg_signals(:,i) = abs(emg_signals(:,i) - mean(emg_signals(:,i)));

% filter EMG data through low pass Butterworth filter 
filt_emg_signals(:,i) = Butterworth_LPF(emg_signals(:,i), 2000, 5, 4);

% normalize EMG to MVC 
norm_emg(:,i) = filt_emg_signals(:,i)/maximum(i) * 100; %percent activation
% Binning data based on type of damping, one if statement for each damping condition and one for each muscle 
if trials(k).TargetDirNum == 4

if trials(k).DampingNumber == 3
    if i == 1
      count_fore_var_4 = count_fore_var_4 + 1;
      fore_var_4(:,count_fore_var_4) = norm_emg(:,i);
    elseif i == 2
      count_bi_var_4 = count_bi_var_4 + 1;
      bi_var_4(:,count_bi_var_4) = norm_emg(:,i);
    elseif i == 3 
      count_front_delt_var_4 = count_front_delt_var_4 + 1;
      front_delt_var_4(:,count_front_delt_var_4) = norm_emg(:,i);
    elseif i == 4 
      count_rear_delt_var_4 = count_rear_delt_var_4 + 1;
      rear_delt_var_4(:,count_rear_delt_var_4) = norm_emg(:,i);
    elseif i == 5
      count_long_tri_var_4 = count_long_tri_var_4 + 1; 
      long_tri_var_4(:,count_long_tri_var_4) = norm_emg(:,i);
    elseif i == 6
      count_lat_tri_var_4 = count_lat_tri_var_4 + 1;
      lat_tri_var_4(:,count_lat_tri_var_4) = norm_emg(:,i);
    end
end

if trials(k).DampingNumber == 2
    if i == 1
      count_fore_neg_4 = count_fore_neg_4 + 1;
      fore_neg_4(:,count_fore_neg_4) = norm_emg(:,i);
    elseif i == 2
      count_bi_neg_4 = count_bi_neg_4 + 1;
      bi_neg_4(:,count_bi_neg_4) = norm_emg(:,i);
    elseif i == 3 
      count_front_delt_neg_4 = count_front_delt_neg_4 + 1;
      front_delt_neg_4(:,count_front_delt_neg_4) = norm_emg(:,i);
    elseif i == 4 
      count_rear_delt_neg_4 = count_rear_delt_neg_4 + 1;
      rear_delt_neg_4(:,count_rear_delt_neg_4) = norm_emg(:,i);
    elseif i == 5
      count_long_tri_neg_4 = count_long_tri_neg_4 + 1; 
      long_tri_neg_4(:,count_long_tri_neg_4) = norm_emg(:,i);
    elseif i == 6
      count_lat_tri_neg_4 = count_lat_tri_neg_4 + 1;
      lat_tri_neg_4(:,count_lat_tri_neg_4) = norm_emg(:,i);
    end
end
if trials(k).DampingNumber == 1
    if i == 1
      count_fore_pos_4 = count_fore_pos_4 + 1;
      fore_pos_4(:,count_fore_pos_4) = norm_emg(:,i);
    elseif i == 2
      count_bi_pos_4 = count_bi_pos_4 + 1;
      bi_pos_4(:,count_bi_pos_4) = norm_emg(:,i);
    elseif i == 3 
      count_front_delt_pos_4 = count_front_delt_pos_4 + 1;
      front_delt_pos_4(:,count_front_delt_pos_4) = norm_emg(:,i);
    elseif i == 4 
      count_rear_delt_pos_4 = count_rear_delt_pos_4 + 1;
      rear_delt_pos_4(:,count_rear_delt_pos_4) = norm_emg(:,i);
    elseif i == 5
      count_long_tri_pos_4 = count_long_tri_pos_4 + 1; 
      long_tri_pos_4(:,count_long_tri_pos_4) = norm_emg(:,i);
    elseif i == 6
      count_lat_tri_pos_4 = count_lat_tri_pos_4 + 1;
      lat_tri_pos_4(:,count_lat_tri_pos_4) = norm_emg(:,i);
    end
end 
end

if trials(k).TargetDirNum == 3

if trials(k).DampingNumber == 3
    if i == 1
      count_fore_var_3 = count_fore_var_3 + 1;
      fore_var_3(:,count_fore_var_3) = norm_emg(:,i);
    elseif i == 2
      count_bi_var_3 = count_bi_var_3 + 1;
      bi_var_3(:,count_bi_var_3) = norm_emg(:,i);
    elseif i == 3 
      count_front_delt_var_3 = count_front_delt_var_3 + 1;
      front_delt_var_3(:,count_front_delt_var_3) = norm_emg(:,i);
    elseif i == 4 
      count_rear_delt_var_3 = count_rear_delt_var_3 + 1;
      rear_delt_var_3(:,count_rear_delt_var_3) = norm_emg(:,i);
    elseif i == 5
      count_long_tri_var_3 = count_long_tri_var_3 + 1; 
      long_tri_var_3(:,count_long_tri_var_3) = norm_emg(:,i);
    elseif i == 6
      count_lat_tri_var_3 = count_lat_tri_var_3 + 1;
      lat_tri_var_3(:,count_lat_tri_var_3) = norm_emg(:,i);
    end
end

if trials(k).DampingNumber == 2
    if i == 1
      count_fore_neg_3 = count_fore_neg_3 + 1;
      fore_neg_3(:,count_fore_neg_3) = norm_emg(:,i);
    elseif i == 2
      count_bi_neg_3 = count_bi_neg_3 + 1;
      bi_neg_3(:,count_bi_neg_3) = norm_emg(:,i);
    elseif i == 3 
      count_front_delt_neg_3 = count_front_delt_neg_3 + 1;
      front_delt_neg_3(:,count_front_delt_neg_3) = norm_emg(:,i);
    elseif i == 4 
      count_rear_delt_neg_3 = count_rear_delt_neg_3 + 1;
      rear_delt_neg_3(:,count_rear_delt_neg_3) = norm_emg(:,i);
    elseif i == 5
      count_long_tri_neg_3 = count_long_tri_neg_3 + 1; 
      long_tri_neg_3(:,count_long_tri_neg_3) = norm_emg(:,i);
    elseif i == 6
      count_lat_tri_neg_3 = count_lat_tri_neg_3 + 1;
      lat_tri_neg_3(:,count_lat_tri_neg_3) = norm_emg(:,i);
    end
end
if trials(k).DampingNumber == 1
    if i == 1
      count_fore_pos_3 = count_fore_pos_3 + 1;
      fore_pos_3(:,count_fore_pos_3) = norm_emg(:,i);
    elseif i == 2
      count_bi_pos_3 =  count_bi_pos_3 + 1;
      bi_pos_3(:,count_bi_pos_3) = norm_emg(:,i);
    elseif i == 3 
      count_front_delt_pos_3 = count_front_delt_pos_3 + 1;
      front_delt_pos_3(:,count_front_delt_pos_3) = norm_emg(:,i);
    elseif i == 4 
      count_rear_delt_pos_3 = count_rear_delt_pos_3 + 1;
      rear_delt_pos_3(:,count_rear_delt_pos_3) = norm_emg(:,i);
    elseif i == 5
      count_long_tri_pos_3 = count_long_tri_pos_3 + 1; 
      long_tri_pos_3(:,count_long_tri_pos_3) = norm_emg(:,i);
    elseif i == 6
      count_lat_tri_pos_3 = count_lat_tri_pos_3 + 1;
      lat_tri_pos_3(:,count_lat_tri_pos_3) = norm_emg(:,i);
    end
end 
end

if trials(k).TargetDirNum == 2
if trials(k).DampingNumber == 3
    if i == 1
      count_fore_var_2 = count_fore_var_2 + 1;
      fore_var_2(:,count_fore_var_2) = norm_emg(:,i);
    elseif i == 2
      count_bi_var_2 = count_bi_var_2 + 1;
      bi_var_2(:,count_bi_var_2) = norm_emg(:,i);
    elseif i == 3 
      count_front_delt_var_2 = count_front_delt_var_2 + 1;
      front_delt_var_2(:,count_front_delt_var_2) = norm_emg(:,i);
    elseif i == 4 
      count_rear_delt_var_2 = count_rear_delt_var_2 + 1;
      rear_delt_var_2(:,count_rear_delt_var_2) = norm_emg(:,i);
    elseif i == 5
      count_long_tri_var_2 = count_long_tri_var_2 + 1; 
      long_tri_var_2(:,count_long_tri_var_2) = norm_emg(:,i);
    elseif i == 6
      count_lat_tri_var_2 = count_lat_tri_var_2 + 1;
      lat_tri_var_2(:,count_lat_tri_var_2) = norm_emg(:,i);
    end
end

if trials(k).DampingNumber == 2
    if i == 1
      count_fore_neg_2 = count_fore_neg_2 + 1;
      fore_neg_2(:,count_fore_neg_2) = norm_emg(:,i);
    elseif i == 2
      count_bi_neg_2 = count_bi_neg_2 + 1;
      bi_neg_2(:,count_bi_neg_2) = norm_emg(:,i);
    elseif i == 3 
      count_front_delt_neg_2 = count_front_delt_neg_2 + 1;
      front_delt_neg_2(:,count_front_delt_neg_2) = norm_emg(:,i);
    elseif i == 4 
      count_rear_delt_neg_2 = count_rear_delt_neg_2 + 1;
      rear_delt_neg_2(:,count_rear_delt_neg_2) = norm_emg(:,i);
    elseif i == 5
      count_long_tri_neg_2 = count_long_tri_neg_2 + 1; 
      long_tri_neg_2(:,count_long_tri_neg_2) = norm_emg(:,i);
    elseif i == 6
      count_lat_tri_neg_2 = count_lat_tri_neg_2 + 1;
      lat_tri_neg_2(:,count_lat_tri_neg_2) = norm_emg(:,i);
    end
end
if trials(k).DampingNumber == 1
    if i == 1
      count_fore_pos_2 = count_fore_pos_2 + 1;
      fore_pos_2(:,count_fore_pos_2) = norm_emg(:,i);
    elseif i == 2
      count_bi_pos_2 =  count_bi_pos_2 + 1;
      bi_pos_2(:,count_bi_pos_2) = norm_emg(:,i);
    elseif i == 3 
      count_front_delt_pos_2 = count_front_delt_pos_2 + 1;
      front_delt_pos_2(:,count_front_delt_pos_2) = norm_emg(:,i);
    elseif i == 4 
      count_rear_delt_pos_2 = count_rear_delt_pos_2 + 1;
      rear_delt_pos_2(:,count_rear_delt_pos_2) = norm_emg(:,i);
    elseif i == 5
      count_long_tri_pos_2 = count_long_tri_pos_2 + 1; 
      long_tri_pos_2(:,count_long_tri_pos_2) = norm_emg(:,i);
    elseif i == 6
      count_lat_tri_pos_2 = count_lat_tri_pos_2 + 1;
      lat_tri_pos_2(:,count_lat_tri_pos_2) = norm_emg(:,i);
    end
end 
end

if trials(k).TargetDirNum == 1
if trials(k).DampingNumber == 3
    if i == 1
      count_fore_var_1 = count_fore_var_1 + 1;
      fore_var_1(:,count_fore_var_1) = norm_emg(:,i);
    elseif i == 2
      count_bi_var_1 = count_bi_var_1 + 1;
      bi_var_1(:,count_bi_var_1) = norm_emg(:,i);
    elseif i == 3 
      count_front_delt_var_1 = count_front_delt_var_1 + 1;
      front_delt_var_1(:,count_front_delt_var_1) = norm_emg(:,i);
    elseif i == 4 
      count_rear_delt_var_1 = count_rear_delt_var_1 + 1;
      rear_delt_var_1(:,count_rear_delt_var_1) = norm_emg(:,i);
    elseif i == 5
      count_long_tri_var_1 = count_long_tri_var_1 + 1; 
      long_tri_var_1(:,count_long_tri_var_1) = norm_emg(:,i);
    elseif i == 6
      count_lat_tri_var_1 = count_lat_tri_var_1 + 1;
      lat_tri_var_1(:,count_lat_tri_var_1) = norm_emg(:,i);
    end
end

if trials(k).DampingNumber == 2
    if i == 1
      count_fore_neg_1 = count_fore_neg_1 + 1;
      fore_neg_1(:,count_fore_neg_1) = norm_emg(:,i);
    elseif i == 2
      count_bi_neg_1 = count_bi_neg_1 + 1;
      bi_neg_1(:,count_bi_neg_1) = norm_emg(:,i);
    elseif i == 3 
      count_front_delt_neg_1 = count_front_delt_neg_1 + 1;
      front_delt_neg_1(:,count_front_delt_neg_1) = norm_emg(:,i);
    elseif i == 4 
      count_rear_delt_neg_1 = count_rear_delt_neg_1 + 1;
      rear_delt_neg_1(:,count_rear_delt_neg_1) = norm_emg(:,i);
    elseif i == 5
      count_long_tri_neg_1 = count_long_tri_neg_1 + 1; 
      long_tri_neg_1(:,count_long_tri_neg_1) = norm_emg(:,i);
    elseif i == 6
      count_lat_tri_neg_1 = count_lat_tri_neg_1 + 1;
      lat_tri_neg_1(:,count_lat_tri_neg_1) = norm_emg(:,i);
    end
end
if trials(k).DampingNumber == 1
    if i == 1
      count_fore_pos_1 = count_fore_pos_1 + 1;
      fore_pos_1(:,count_fore_pos_1) = norm_emg(:,i);
    elseif i == 2
      count_bi_pos_1 =  count_bi_pos_1 + 1;
      bi_pos_1(:,count_bi_pos_1) = norm_emg(:,i);
    elseif i == 3 
      count_front_delt_pos_1 = count_front_delt_pos_1 + 1;
      front_delt_pos_1(:,count_front_delt_pos_1) = norm_emg(:,i);
    elseif i == 4 
      count_rear_delt_pos_1 = count_rear_delt_pos_1 + 1;
      rear_delt_pos_1(:,count_rear_delt_pos_1) = norm_emg(:,i);
    elseif i == 5
      count_long_tri_pos_1 = count_long_tri_pos_1 + 1; 
      long_tri_pos_1(:,count_long_tri_pos_1) = norm_emg(:,i);
    elseif i == 6
      count_lat_tri_pos_1 = count_lat_tri_pos_1 + 1;
      lat_tri_pos_1(:,count_lat_tri_pos_1) = norm_emg(:,i);
    end
end 
end

end 
end
end
end


%% 7) Mean Calculations
%Calculating average of each vector, needed to transpose in-order for mean to work correctly
%Averages for upward direction 

% fore_avg_var_4 = mean(fore_var_4');
% fore_avg_pos_4 = mean(fore_pos_4');
% fore_avg_neg_4 = mean(fore_neg_4');
% 
% bi_avg_var_4 = mean(bi_var_4');
% bi_avg_pos_4 = mean(bi_pos_4');
% bi_avg_neg_4 = mean(bi_neg_4');
% 
% front_delt_avg_var_4 = mean(front_delt_var_4');
% front_delt_avg_pos_4 = mean(front_delt_pos_4');
% front_delt_avg_neg_4 = mean(front_delt_neg_4');
% 
% rear_delt_avg_var_4 = mean(rear_delt_var_4');
% rear_delt_avg_pos_4 = mean(rear_delt_pos_4');
% rear_delt_avg_neg_4 = mean(rear_delt_neg_4');
% 
% long_tri_avg_var_4 = mean(long_tri_var_4');
% long_tri_avg_pos_4 = mean(long_tri_pos_4');
% long_tri_avg_neg_4 = mean(long_tri_neg_4');
% 
% lat_tri_avg_var_4 = mean(lat_tri_var_4');
% lat_tri_avg_pos_4 = mean(lat_tri_pos_4');
% lat_tri_avg_neg_4 = mean(lat_tri_neg_4');
% 
% % Average activation for upwards 
% mean_4_var = mean((fore_avg_var_4+bi_avg_var_4+front_delt_avg_var_4+rear_delt_avg_var_4+long_tri_avg_var_4+lat_tri_avg_var_4)/6);
% mean_4_pos = mean((fore_avg_pos_4+bi_avg_pos_4+front_delt_avg_pos_4+rear_delt_avg_pos_4+long_tri_avg_pos_4+lat_tri_avg_pos_4)/6);
% mean_4_neg = mean((fore_avg_neg_4+bi_avg_neg_4+bi_avg_neg_4+front_delt_avg_neg_4+rear_delt_avg_neg_4+long_tri_avg_neg_4+lat_tri_avg_neg_4)/6);
% 
% % Individual Muscle Activation 
% mfore_avg_var_4 = mean(fore_avg_var_4);
% mfore_avg_pos_4 = mean(fore_avg_pos_4);
% mfore_avg_neg_4 = mean(fore_avg_neg_4);
% 
% mbi_avg_var_4 = mean(bi_avg_var_4);
% mbi_avg_pos_4 = mean(bi_avg_pos_4);
% mbi_avg_neg_4 = mean(bi_avg_neg_4);
% 
% mfront_delt_avg_pos_4 = mean(front_delt_avg_pos_4);
% mfront_delt_avg_var_4 = mean(front_delt_avg_var_4);
% mfront_delt_avg_neg_4 = mean(front_delt_avg_neg_4);
% 
% mrear_delt_avg_pos_4 = mean(rear_delt_avg_pos_4);
% mrear_delt_avg_var_4 = mean(rear_delt_avg_var_4);
% mrear_delt_avg_neg_4 = mean(rear_delt_avg_neg_4);
% 
% mlong_tri_avg_pos_4 = mean(long_tri_avg_pos_4);
% mlong_tri_avg_var_4 = mean(long_tri_avg_var_4);
% mlong_tri_avg_neg_4 = mean(long_tri_avg_neg_4);
% 
% mlat_tri_avg_pos_4 = mean(lat_tri_avg_pos_4);
% mlat_tri_avg_var_4 = mean(lat_tri_avg_var_4);
% mlat_tri_avg_neg_4 = mean(lat_tri_avg_neg_4);
% 
% % Calculating averages for downward direction 
% fore_avg_var_3 = mean(fore_var_3');
% fore_avg_pos_3 = mean(fore_pos_3');
% fore_avg_neg_3 = mean(fore_neg_3');
% 
% bi_avg_var_3 = mean(bi_var_3');
% bi_avg_pos_3 = mean(bi_pos_3');
% bi_avg_neg_3 = mean(bi_neg_3');
% 
% front_delt_avg_var_3 = mean(front_delt_var_3');
% front_delt_avg_pos_3 = mean(front_delt_pos_3');
% front_delt_avg_neg_3 = mean(front_delt_neg_3');
% 
% rear_delt_avg_var_3 = mean(rear_delt_var_3');
% rear_delt_avg_pos_3 = mean(rear_delt_pos_3');
% rear_delt_avg_neg_3 = mean(rear_delt_neg_3');
% 
% long_tri_avg_var_3 = mean(long_tri_var_3');
% long_tri_avg_pos_3 = mean(long_tri_pos_3');
% long_tri_avg_neg_3 = mean(long_tri_neg_3');
% 
% lat_tri_avg_var_3 = mean(lat_tri_var_3');
% lat_tri_avg_pos_3 = mean(lat_tri_pos_3');
% lat_tri_avg_neg_3 = mean(lat_tri_neg_3');
% 
% mean_3_var = mean((fore_avg_var_3+bi_avg_neg_3+front_delt_avg_neg_3+rear_delt_avg_neg_3+long_tri_avg_neg_3+lat_tri_avg_neg_3)/6);
% mean_3_pos = mean((fore_avg_pos_3+bi_avg_pos_3+front_delt_avg_pos_3+rear_delt_avg_pos_3+long_tri_avg_pos_3+lat_tri_avg_pos_3)/6);
% mean_3_neg = mean((fore_avg_neg_3+bi_avg_neg_3+front_delt_avg_neg_3+rear_delt_avg_neg_3+long_tri_avg_neg_3+lat_tri_avg_neg_3)/6);
% 
% % Individual Muscle Activation Downwards
% mfore_avg_var_3 = mean(fore_avg_var_3);
% mfore_avg_pos_3 = mean(fore_avg_pos_3);
% mfore_avg_neg_3 = mean(fore_avg_neg_3);
% 
% mbi_avg_var_3 = mean(bi_avg_var_3);
% mbi_avg_pos_3 = mean(bi_avg_pos_3);
% mbi_avg_neg_3 = mean(bi_avg_neg_3);
% 
% mfront_delt_avg_pos_3 = mean(front_delt_avg_pos_3);
% mfront_delt_avg_var_3 = mean(front_delt_avg_var_3);
% mfront_delt_avg_neg_3 = mean(front_delt_avg_neg_3);
% 
% mrear_delt_avg_pos_3 = mean(rear_delt_avg_pos_3);
% mrear_delt_avg_var_3 = mean(rear_delt_avg_var_3);
% mrear_delt_avg_neg_3 = mean(rear_delt_avg_neg_3);
% 
% mlong_tri_avg_pos_3 = mean(long_tri_avg_pos_3);
% mlong_tri_avg_var_3 = mean(long_tri_avg_var_3);
% mlong_tri_avg_neg_3 = mean(long_tri_avg_neg_3);
% 
% mlat_tri_avg_pos_3 = mean(lat_tri_avg_pos_3);
% mlat_tri_avg_var_3 = mean(lat_tri_avg_var_3);
% mlat_tri_avg_neg_3 = mean(lat_tri_avg_neg_3);
% % 
% Calculating averages for right direction 
fore_avg_var_2 = mean(fore_var_2');
fore_avg_pos_2 = mean(fore_pos_2');
fore_avg_neg_2 = mean(fore_neg_2');

bi_avg_var_2 = mean(bi_var_2');
bi_avg_pos_2 = mean(bi_pos_2');
bi_avg_neg_2 = mean(bi_neg_2');

front_delt_avg_var_2 = mean(front_delt_var_2');
front_delt_avg_pos_2 = mean(front_delt_pos_2');
front_delt_avg_neg_2 = mean(front_delt_neg_2');

rear_delt_avg_var_2 = mean(rear_delt_var_2');
rear_delt_avg_pos_2 = mean(rear_delt_pos_2');
rear_delt_avg_neg_2 = mean(rear_delt_neg_2');

long_tri_avg_var_2 = mean(long_tri_var_2');
long_tri_avg_pos_2 = mean(long_tri_pos_2');
long_tri_avg_neg_2 = mean(long_tri_neg_2');

lat_tri_avg_var_2 = mean(lat_tri_var_2');
lat_tri_avg_pos_2 = mean(lat_tri_pos_2');
lat_tri_avg_neg_2 = mean(lat_tri_neg_2');

% mean_2_var = mean((fore_avg_var_2+bi_avg_var_2+front_delt_avg_var_2+rear_delt_avg_var_2+long_tri_avg_var_2+lat_tri_avg_var_2)/6);
% mean_2_pos = mean((fore_avg_pos_2+bi_avg_pos_2+front_delt_avg_pos_2+rear_delt_avg_pos_2+long_tri_avg_pos_2+lat_tri_avg_pos_2)/6);
% mean_2_neg = mean((fore_avg_neg_2+bi_avg_neg_2+front_delt_avg_neg_2+rear_delt_avg_neg_2+long_tri_avg_neg_2+lat_tri_avg_neg_2)/6);

% Individual Muscle Activation Right
mfore_avg_var_2 = mean(fore_avg_var_2);
mfore_avg_pos_2 = mean(fore_avg_pos_2);
mfore_avg_neg_2 = mean(fore_avg_neg_2);

mbi_avg_var_2 = mean(bi_avg_var_2);
mbi_avg_pos_2 = mean(bi_avg_pos_2);
mbi_avg_neg_2 = mean(bi_avg_neg_2);

mfront_delt_avg_pos_2 = mean(front_delt_avg_pos_2);
mfront_delt_avg_var_2 = mean(front_delt_avg_var_2);
mfront_delt_avg_neg_2 = mean(front_delt_avg_neg_2);

mrear_delt_avg_pos_2 = mean(rear_delt_avg_pos_2);
mrear_delt_avg_var_2 = mean(rear_delt_avg_var_2);
mrear_delt_avg_neg_2 = mean(rear_delt_avg_neg_2);

mlong_tri_avg_pos_2 = mean(long_tri_avg_pos_2);
mlong_tri_avg_var_2 = mean(long_tri_avg_var_2);
mlong_tri_avg_neg_2 = mean(long_tri_avg_neg_2);

mlat_tri_avg_pos_2 = mean(lat_tri_avg_pos_2);
mlat_tri_avg_var_2 = mean(lat_tri_avg_var_2);
mlat_tri_avg_neg_2 = mean(lat_tri_avg_neg_2);


% Calculating averages for left direction 
fore_avg_var_1 = mean(fore_var_1');
fore_avg_pos_1 = mean(fore_pos_1');
fore_avg_neg_1 = mean(fore_neg_1');

bi_avg_var_1 = mean(bi_var_1');
bi_avg_pos_1 = mean(bi_pos_1');
bi_avg_neg_1 = mean(bi_neg_1');

front_delt_avg_var_1 = mean(front_delt_var_1');
front_delt_avg_pos_1 = mean(front_delt_pos_1');
front_delt_avg_neg_1 = mean(front_delt_neg_1');

rear_delt_avg_var_1 = mean(rear_delt_var_1');
rear_delt_avg_pos_1 = mean(rear_delt_pos_1');
rear_delt_avg_neg_1 = mean(rear_delt_neg_1');

long_tri_avg_var_1 = mean(long_tri_var_1');
long_tri_avg_pos_1 = mean(long_tri_pos_1');
long_tri_avg_neg_1 = mean(long_tri_neg_1');

lat_tri_avg_var_1 = mean(lat_tri_var_1');
lat_tri_avg_pos_1 = mean(lat_tri_pos_1');
lat_tri_avg_neg_1 = mean(lat_tri_neg_1');

% mean_1_var = mean((fore_avg_var_1+bi_avg_var_1+front_delt_avg_var_1+rear_delt_avg_var_1+long_tri_avg_var_1+lat_tri_avg_var_1)/6);
% mean_1_pos = mean((fore_avg_pos_1+bi_avg_pos_1+front_delt_avg_pos_1+rear_delt_avg_pos_1+long_tri_avg_pos_1+lat_tri_avg_pos_1)/6);
% mean_1_neg = mean((fore_avg_neg_1+bi_avg_neg_1+front_delt_avg_neg_1+rear_delt_avg_neg_1+long_tri_avg_neg_1+lat_tri_avg_neg_1)/6);

% Individual Muscle Activation Left
mfore_avg_var_1 = mean(fore_avg_var_1);
mfore_avg_pos_1 = mean(fore_avg_pos_1);
mfore_avg_neg_1 = mean(fore_avg_neg_1);

mbi_avg_var_1 = mean(bi_avg_var_1);
mbi_avg_pos_1 = mean(bi_avg_pos_1);
mbi_avg_neg_1 = mean(bi_avg_neg_1);

mfront_delt_avg_pos_1 = mean(front_delt_avg_pos_1);
mfront_delt_avg_var_1 = mean(front_delt_avg_var_1);
mfront_delt_avg_neg_1 = mean(front_delt_avg_neg_1);

mrear_delt_avg_pos_1 = mean(rear_delt_avg_pos_1);
mrear_delt_avg_var_1 = mean(rear_delt_avg_var_1);
mrear_delt_avg_neg_1 = mean(rear_delt_avg_neg_1);

mlong_tri_avg_pos_1 = mean(long_tri_avg_pos_1);
mlong_tri_avg_var_1 = mean(long_tri_avg_var_1);
mlong_tri_avg_neg_1 = mean(long_tri_avg_neg_1);

mlat_tri_avg_pos_1 = mean(lat_tri_avg_pos_1);
mlat_tri_avg_var_1 = mean(lat_tri_avg_var_1);
mlat_tri_avg_neg_1 = mean(lat_tri_avg_neg_1);

%% 8) Plotting Percent Activation

% Upwards Dir.
figure;
subplot (3,2,1)
plot(x_axis,fore_avg_pos_4,x_axis,fore_avg_neg_4,x_axis,fore_avg_var_4); title("Muscle Activation of Forearm for Upward Direction"); xlabel("Time (s)"); ylabel("Percent Activation");legend('Positive','Negative','Variable');
subplot (3,2,2)
plot(x_axis,bi_avg_pos_4,x_axis,bi_avg_neg_4,x_axis,bi_avg_var_4); title("Muscle Activation of Bicep for Upward Direction"); xlabel("Time (s)"); ylabel("Percent Activation");legend('Positive','Negative','Variable');
subplot (3,2,3)
plot(x_axis,front_delt_avg_pos_4,x_axis,front_delt_avg_neg_4,x_axis,front_delt_avg_var_4); title("Muscle Activation of Front Delt. for Upward Direction"); xlabel("Time (s)"); ylabel("Percent Activation");legend('Positive','Negative','Variable');
subplot (3,2,4)
plot(x_axis,rear_delt_avg_pos_4,x_axis,rear_delt_avg_neg_4,x_axis,rear_delt_avg_var_4); title("Muscle Activation of Rear Delt. for Upward Direction"); xlabel("Time (s)"); ylabel("Percent Activation");legend('Positive','Negative','Variable');
subplot (3,2,5) 
plot(x_axis,long_tri_avg_pos_4,x_axis,long_tri_avg_neg_4,x_axis,long_tri_avg_var_4); title("Muscle Activation of Long. Tri. for Upward Direction"); xlabel("Time (s)"); ylabel("Percent Activation");legend('Positive','Negative','Variable');
subplot (3,2,6) 
plot(x_axis,lat_tri_avg_pos_4,x_axis,lat_tri_avg_neg_4,x_axis,lat_tri_avg_var_4); title("Muscle Activation of Lat. Tri. for Upward Direction"); xlabel("Time (s)"); ylabel("Percent Activation");legend('Positive','Negative','Variable');

% Downwards Dir.
figure;
subplot (3,2,1)
plot(x_axis,fore_avg_pos_3,x_axis,fore_avg_neg_3,x_axis,fore_avg_var_3); title("Muscle Activation of Forearm for Downward Direction"); xlabel("Time (s)"); ylabel("Percent Activation");legend('Positive','Negative','Variable');
subplot (3,2,2)
plot(x_axis,bi_avg_pos_3,x_axis,bi_avg_neg_3,x_axis,bi_avg_var_3); title("Muscle Activation of Bicep for Downward Direction"); xlabel("Time (s)"); ylabel("Percent Activation");legend('Positive','Negative','Variable');
subplot (3,2,3)
plot(x_axis,front_delt_avg_pos_3,x_axis,front_delt_avg_neg_3,x_axis,front_delt_avg_var_3); title("Muscle Activation of Front Delt. for Downward Direction"); xlabel("Time (s)"); ylabel("Percent Activation");legend('Positive','Negative','Variable');
subplot (3,2,4)
plot(x_axis,rear_delt_avg_pos_3,x_axis,rear_delt_avg_neg_3,x_axis,rear_delt_avg_var_3); title("Muscle Activation of Rear Delt. for Downward Direction"); xlabel("Time (s)"); ylabel("Percent Activation");legend('Positive','Negative','Variable');
subplot (3,2,5)
plot(x_axis,long_tri_avg_pos_3,x_axis,long_tri_avg_neg_3,x_axis,long_tri_avg_var_3); title("Muscle Activation of Long. Tri. for Downward Direction"); xlabel("Time (s)"); ylabel("Percent Activation");legend('Positive','Negative','Variable');
subplot (3,2,6)
plot(x_axis,lat_tri_avg_pos_3,x_axis,lat_tri_avg_neg_3,x_axis,lat_tri_avg_var_3); title("Muscle Activation of Lat. Tri. for Downward Direction"); xlabel("Time (s)"); ylabel("Percent Activation");legend('Positive','Negative','Variable');
% 
% % Right Direction
% figure;
% subplot (3,2,1)
% plot(x_axis,fore_avg_pos_2,x_axis,fore_avg_neg_2,x_axis,fore_avg_var_2); title("Muscle Activation of Forearm for Right Direction"); xlabel("Time (s)"); ylabel("Percent Activation");legend('Positive','Negative','Variable');
% subplot (3,2,2)
% plot(x_axis,bi_avg_pos_2,x_axis,bi_avg_neg_2,x_axis,bi_avg_var_2); title("Muscle Activation of Bicep for Right Direction"); xlabel("Time (s)"); ylabel("Percent Activation");legend('Positive','Negative','Variable');
% subplot (3,2,3)
% plot(x_axis,front_delt_avg_pos_2,x_axis,front_delt_avg_neg_2,x_axis,front_delt_avg_var_2); title("Muscle Activation of Front Delt. for Right Direction"); xlabel("Time (s)"); ylabel("Percent Activation");legend('Positive','Negative','Variable');
% subplot (3,2,4)
% plot(x_axis,rear_delt_avg_pos_2,x_axis,rear_delt_avg_neg_2,x_axis,rear_delt_avg_var_2); title("Muscle Activation of Rear Delt. for Right Direction"); xlabel("Time (s)"); ylabel("Percent Activation");legend('Positive','Negative','Variable');
% subplot (3,2,5)
% plot(x_axis,long_tri_avg_pos_2,x_axis,long_tri_avg_neg_2,x_axis,long_tri_avg_var_2); title("Muscle Activation of Long. Tri. for Right Direction"); xlabel("Time (s)"); ylabel("Percent Activation");legend('Positive','Negative','Variable');
% subplot (3,2,6)
% plot(x_axis,lat_tri_avg_pos_2,x_axis,lat_tri_avg_neg_2,x_axis,lat_tri_avg_var_2); title("Muscle Activation of Lat. Tri. for Right Direction"); xlabel("Time (s)"); ylabel("Percent Activation");legend('Positive','Negative','Variable');
% 
% % Left Direction 
% figure;
% subplot (3,2,1)
% plot(x_axis,fore_avg_pos_1,x_axis,fore_avg_neg_1,x_axis,fore_avg_var_1); title("Muscle Activation of Forearm for Left Direction"); xlabel("Time (s)"); ylabel("Percent Activation");legend('Positive','Negative','Variable');
% subplot (3,2,2)
% plot(x_axis,bi_avg_pos_1,x_axis,bi_avg_neg_1,x_axis,bi_avg_var_1); title("Muscle Activation of Bicep for Left Direction"); xlabel("Time (s)"); ylabel("Percent Activation");legend('Positive','Negative','Variable');
% subplot (3,2,3)
% plot(x_axis,front_delt_avg_pos_1,x_axis,front_delt_avg_neg_1,x_axis,front_delt_avg_var_1); title("Muscle Activation of Front Delt. for Left Direction"); xlabel("Time (s)"); ylabel("Percent Activation");legend('Positive','Negative','Variable');
% subplot (3,2,4)
% plot(x_axis,rear_delt_avg_pos_1,x_axis,rear_delt_avg_neg_1,x_axis,rear_delt_avg_var_1); title("Muscle Activation of Rear Delt. for Left Direction"); xlabel("Time (s)"); ylabel("Percent Activation");legend('Positive','Negative','Variable');
% subplot (3,2,5)
% plot(x_axis,long_tri_avg_pos_1,x_axis,long_tri_avg_neg_1,x_axis,long_tri_avg_var_1); title("Muscle Activation of Long. Tri. for Left Direction"); xlabel("Time (s)"); ylabel("Percent Activation");legend('Positive','Negative','Variable');
% subplot (3,2,6)
% plot(x_axis,lat_tri_avg_pos_1,x_axis,lat_tri_avg_neg_1,x_axis,lat_tri_avg_var_1); title("Muscle Activation of Lat. Tri. for Left Direction"); xlabel("Time (s)"); ylabel("Percent Activation");legend('Positive','Negative','Variable');

%% Movement Phase Muscle Activation Per Direction 
% Upwards 
figure; 
subplot (1,2,1);
y = [ mfore_avg_pos_4 mfore_avg_var_4 mfore_avg_neg_4; mbi_avg_pos_4 mbi_avg_var_4 mbi_avg_neg_4; mfront_delt_avg_pos_4 mfront_delt_avg_var_4 mfront_delt_avg_neg_4; mrear_delt_avg_pos_4 mrear_delt_avg_var_4 mrear_delt_avg_neg_4; mlong_tri_avg_pos_4 mlong_tri_avg_var_4 mlong_tri_avg_neg_4; mlat_tri_avg_pos_4 mlat_tri_avg_var_4 mlat_tri_avg_neg_4 ];
x = categorical({'Forearm','Bicep','Front Delt','Rear Delt', 'Longitudinal Tricep', 'Lateral Tricep'});
x = reordercats(x,{'Forearm','Bicep','Front Delt','Rear Delt', 'Longitudinal Tricep', 'Lateral Tricep'});

bar(x,y); legend('Positive','Variable','Negative'); title("Average Muscle Activation for Movement Period in Upwards Dir."); ylabel("Percent Activation");xlabel("Muscle");

% Downwards 
subplot(1,2,2);
y = [ mfore_avg_pos_3 mfore_avg_var_3 mfore_avg_neg_3; mbi_avg_pos_3 mbi_avg_var_3 mbi_avg_neg_3; mfront_delt_avg_pos_3 mfront_delt_avg_var_3 mfront_delt_avg_neg_3; mrear_delt_avg_pos_3 mrear_delt_avg_var_3 mrear_delt_avg_neg_3; mlong_tri_avg_pos_3 mlong_tri_avg_var_3 mlong_tri_avg_neg_3; mlat_tri_avg_pos_3 mlat_tri_avg_var_3 mlat_tri_avg_neg_3 ];
x = categorical({'Forearm','Bicep','Front Delt','Rear Delt', 'Longitudinal Tricep', 'Lateral Tricep'});
x = reordercats(x,{'Forearm','Bicep','Front Delt','Rear Delt', 'Longitudinal Tricep', 'Lateral Tricep'});

bar(x,y); legend('Positive','Variable','Negative'); title("Average Muscle Activation for Movement Period in Downwards Dir."); ylabel("Percent Activation");xlabel("Muscle");

% % Right 
% subplot(2,1,1);
% y = [ mfore_avg_pos_2 mfore_avg_var_2 mfore_avg_neg_2; mbi_avg_pos_2 mbi_avg_var_2 mbi_avg_neg_2; mfront_delt_avg_pos_2 mfront_delt_avg_var_2 mfront_delt_avg_neg_2; mrear_delt_avg_pos_2 mrear_delt_avg_var_2 mrear_delt_avg_neg_2; mlong_tri_avg_pos_2 mlong_tri_avg_var_2 mlong_tri_avg_neg_2; mlat_tri_avg_pos_2 mlat_tri_avg_var_2 mlat_tri_avg_neg_2 ];
% x = categorical({'Forearm','Bicep','Front Delt','Rear Delt', 'Longitudinal Tricep', 'Lateral Tricep'});
% x = reordercats(x,{'Forearm','Bicep','Front Delt','Rear Delt', 'Longitudinal Tricep', 'Lateral Tricep'});
% 
% bar(x,y); legend('Positive','Variable','Negative'); title("Average Muscle Activation for Movement Period in Right Dir."); ylabel("Percent Activation");xlabel("Muscle");
% 
% % Left 
% subplot(2,1,2);
% y = [ mfore_avg_pos_1 mfore_avg_var_1 mfore_avg_neg_1; mbi_avg_pos_1 mbi_avg_var_1 mbi_avg_neg_1; mfront_delt_avg_pos_1 mfront_delt_avg_var_1 mfront_delt_avg_neg_1; mrear_delt_avg_pos_1 mrear_delt_avg_var_1 mrear_delt_avg_neg_1; mlong_tri_avg_pos_1 mlong_tri_avg_var_1 mlong_tri_avg_neg_1; mlat_tri_avg_pos_1 mlat_tri_avg_var_1 mlat_tri_avg_neg_1 ];
% x = categorical({'Forearm','Bicep','Front Delt','Rear Delt', 'Longitudinal Tricep', 'Lateral Tricep'});
% x = reordercats(x,{'Forearm','Bicep','Front Delt','Rear Delt', 'Longitudinal Tricep', 'Lateral Tricep'});
% 
% bar(x,y); legend('Positive','Variable','Negative'); title("Average Muscle Activation for Movement Period in Left Dir."); ylabel("Percent Activation");xlabel("Muscle");
% 

%%  Mean Muscle Activation Across Directions for Movement Phase 
figure; 
y = [ mean_1_pos mean_1_var mean_1_neg;  mean_2_pos mean_2_var mean_2_neg ; mean_3_pos mean_3_var mean_3_neg; mean_4_pos mean_4_var mean_4_neg];
x = categorical({'Left','Right','Upwards','Downwards'});
x = reordercats(x,{'Left','Right','Upwards','Downwards'});
bar(x,y); legend('Positive','Variable','Negative'); title("Average Muscle Activation Across Directions for Movement Period"); ylabel("Percent Activation");xlabel("Direction");

%% Stability Phase Muscle Activation Per Direction 
% Upwards 
figure; 
subplot (2,2,1);
y = [ mfore_avg_pos_4 mfore_avg_var_4 mfore_avg_neg_4; mbi_avg_pos_4 mbi_avg_var_4 mbi_avg_neg_4; mfront_delt_avg_pos_4 mfront_delt_avg_var_4 mfront_delt_avg_neg_4; mrear_delt_avg_pos_4 mrear_delt_avg_var_4 mrear_delt_avg_neg_4; mlong_tri_avg_pos_4 mlong_tri_avg_var_4 mlong_tri_avg_neg_4; mlat_tri_avg_pos_4 mlat_tri_avg_var_4 mlat_tri_avg_neg_4 ];
x = categorical({'Forearm','Bicep','Front Delt','Rear Delt', 'Longitudinal Tricep', 'Lateral Tricep'});
x = reordercats(x,{'Forearm','Bicep','Front Delt','Rear Delt', 'Longitudinal Tricep', 'Lateral Tricep'});

bar(x,y); legend('Positive','Variable','Negative'); title("Average Muscle Activation for Stability Period in Upwards Dir."); ylabel("Percent Activation");xlabel("Muscle");

% Downwards 
subplot(2,2,2);
y = [ mfore_avg_pos_3 mfore_avg_var_3 mfore_avg_neg_3; mbi_avg_pos_3 mbi_avg_var_3 mbi_avg_neg_3; mfront_delt_avg_pos_3 mfront_delt_avg_var_3 mfront_delt_avg_neg_3; mrear_delt_avg_pos_3 mrear_delt_avg_var_3 mrear_delt_avg_neg_3; mlong_tri_avg_pos_3 mlong_tri_avg_var_3 mlong_tri_avg_neg_3; mlat_tri_avg_pos_3 mlat_tri_avg_var_3 mlat_tri_avg_neg_3 ];
x = categorical({'Forearm','Bicep','Front Delt','Rear Delt', 'Longitudinal Tricep', 'Lateral Tricep'});
x = reordercats(x,{'Forearm','Bicep','Front Delt','Rear Delt', 'Longitudinal Tricep', 'Lateral Tricep'});

bar(x,y); legend('Positive','Variable','Negative'); title("Average Muscle Activation for Stability Period in Downwards Dir."); ylabel("Percent Activation");xlabel("Muscle");

% Right 
subplot(2,2,3);
y = [ mfore_avg_pos_2 mfore_avg_var_2 mfore_avg_neg_2; mbi_avg_pos_2 mbi_avg_var_2 mbi_avg_neg_2; mfront_delt_avg_pos_2 mfront_delt_avg_var_2 mfront_delt_avg_neg_2; mrear_delt_avg_pos_2 mrear_delt_avg_var_2 mrear_delt_avg_neg_2; mlong_tri_avg_pos_2 mlong_tri_avg_var_2 mlong_tri_avg_neg_2; mlat_tri_avg_pos_2 mlat_tri_avg_var_2 mlat_tri_avg_neg_2 ];
x = categorical({'Forearm','Bicep','Front Delt','Rear Delt', 'Longitudinal Tricep', 'Lateral Tricep'});
x = reordercats(x,{'Forearm','Bicep','Front Delt','Rear Delt', 'Longitudinal Tricep', 'Lateral Tricep'});

bar(x,y); legend('Positive','Variable','Negative'); title("Average Muscle Activation for Stability Period in Right Dir."); ylabel("Percent Activation");xlabel("Muscle");

% Left 
subplot(2,2,4);
y = [ mfore_avg_pos_1 mfore_avg_var_1 mfore_avg_neg_1; mbi_avg_pos_1 mbi_avg_var_1 mbi_avg_neg_1; mfront_delt_avg_pos_1 mfront_delt_avg_var_1 mfront_delt_avg_neg_1; mrear_delt_avg_pos_1 mrear_delt_avg_var_1 mrear_delt_avg_neg_1; mlong_tri_avg_pos_1 mlong_tri_avg_var_1 mlong_tri_avg_neg_1; mlat_tri_avg_pos_1 mlat_tri_avg_var_1 mlat_tri_avg_neg_1 ];
x = categorical({'Forearm','Bicep','Front Delt','Rear Delt', 'Longitudinal Tricep', 'Lateral Tricep'});
x = reordercats(x,{'Forearm','Bicep','Front Delt','Rear Delt', 'Longitudinal Tricep', 'Lateral Tricep'});

bar(x,y); legend('Positive','Variable','Negative'); title("Average Muscle Activation for Stability Period in Left Dir."); ylabel("Percent Activation");xlabel("Muscle");

%%  Mean Muscle Activation Across Directions for Stability 
figure; 
y = [ mean_1_pos mean_1_var mean_1_neg;  mean_2_pos mean_2_var mean_2_neg ; mean_3_pos mean_3_var mean_3_neg; mean_4_pos mean_4_var mean_4_neg];
x = categorical({'Left','Right','Upwards','Downwards'});
x = reordercats(x,{'Left','Right','Upwards','Downwards'});
bar(x,y); legend('Positive','Variable','Negative'); title("Average Muscle Activation Across Directions for Stability Period"); ylabel("Percent Activation");xlabel("Direction");

%% 2) Plotting MVC
figure;
for i=1:6
% MVC Comparison of Filtered and Unfiltered Data 
subplot(3,2,i)
muscles = ["Forearm","Bicep","Front Delt","Rear Delt","Longitudinal Tricep","Lateral Tricep"];
title_mvc_1 = sprintf('MVC Signal Comparison for %s',muscles(i));
plot(x_axis_mvc, rect_mvc(:,i));xlabel('Time (s)');ylabel('Signal Amplitude (mV)'); 
hold on;
plot(x_axis_mvc, filt_mvc_signals(:,i),'r','linewidth',2); title(title_mvc_1);xlabel('Time (s)');ylabel('Signal Amplitude (V)'); legend('Unfiltered','Filtered');
end

%% 
x = [ mfront_delt_avg_pos_2 mfront_delt_avg_var_2 mfront_delt_avg_neg_2];
bar(x);
%% Plot for FURI Poster 
% For Right 
% Agonists are: Forearm, Rear Delt, LongT, LatT
figure;
subplot(1,2,1);
y = [ mfore_avg_pos_2 mfore_avg_var_2 mfore_avg_neg_2; mrear_delt_avg_pos_2 mrear_delt_avg_var_2 mrear_delt_avg_neg_2; mlong_tri_avg_pos_2 mlong_tri_avg_var_2 mlong_tri_avg_neg_2; mlat_tri_avg_pos_2 mlat_tri_avg_var_2 mlat_tri_avg_neg_2 ];
x = categorical({'Forearm','Rear Delt', 'Long. Tricep', 'Lateral Tricep'});
x = reordercats(x,{'Forearm','Rear Delt', 'Long. Tricep', 'Lateral Tricep'});

bar(x,y); legend('Positive','Variable','Negative'); title("Average Muscle Activation for Movement Phase in Right Dir."); ylabel("Percent Activation");xlabel("Muscle");

% For Left 
% Agonists are: Forearm, Bicep, LongT, LatT 
subplot(1,2,2);
y = [ mfore_avg_pos_1 mfore_avg_var_1 mfore_avg_neg_1; mbi_avg_pos_1 mbi_avg_var_1 mbi_avg_neg_1; mfront_delt_avg_pos_1 mfront_delt_avg_var_1 mfront_delt_avg_neg_1; mlong_tri_avg_pos_1 mlong_tri_avg_var_1 mlong_tri_avg_neg_1]; 
x = categorical({'Forearm','Bicep','Front Deltoid' 'Long. Tricep'});
x = reordercats(x,{'Forearm','Bicep','Front Deltoid', 'Long. Tricep'});

bar(x,y); legend('Positive','Variable','Negative'); title("Average Muscle Activation for Movement Phase in Left Dir."); ylabel("Percent Activation");xlabel("Muscle");
 
%% Stability Phase 
% For Right 
% Agonists are: Forearm, Rear Delt, LongT, LatT
figure;
subplot(1,2,1);
y = [ mfore_avg_pos_2 mfore_avg_var_2 mfore_avg_neg_2; mrear_delt_avg_pos_2 mrear_delt_avg_var_2 mrear_delt_avg_neg_2; mlong_tri_avg_pos_2 mlong_tri_avg_var_2 mlong_tri_avg_neg_2; mlat_tri_avg_pos_2 mlat_tri_avg_var_2 mlat_tri_avg_neg_2 ];
x = categorical({'Forearm','Rear Delt', 'Long. Tricep', 'Lateral Tricep'});
x = reordercats(x,{'Forearm','Rear Delt', 'Long. Tricep', 'Lateral Tricep'});

bar(x,y); legend('Positive','Variable','Negative'); title("Average Muscle Activation for Stability Phase in Right Dir."); ylabel("Percent Activation");xlabel("Muscle");

% For Left 
% Agonists are: Forearm, Bicep, LongT, LatT 
subplot(1,2,2);
y = [ mfore_avg_pos_1 mfore_avg_var_1 mfore_avg_neg_1; mbi_avg_pos_1 mbi_avg_var_1 mbi_avg_neg_1; mfront_delt_avg_pos_1 mfront_delt_avg_var_1 mfront_delt_avg_neg_1; mlong_tri_avg_pos_1 mlong_tri_avg_var_1 mlong_tri_avg_neg_1]; 
x = categorical({'Forearm','Bicep','Front Deltoid' 'Long. Tricep'});
x = reordercats(x,{'Forearm','Bicep','Front Deltoid', 'Long. Tricep'});

bar(x,y); legend('Positive','Variable','Negative'); title("Average Muscle Activation for Stability Phase in Left Dir."); ylabel("Percent Activation");xlabel("Muscle");

%% FURI Plot 2 
% Upwards 
% Agonists are triceps and front delt
figure; 
subplot (1,2,1);
y = [mfront_delt_avg_pos_4 mfront_delt_avg_var_4 mfront_delt_avg_neg_4; mlong_tri_avg_pos_4 mlong_tri_avg_var_4 mlong_tri_avg_neg_4; mlat_tri_avg_pos_4 mlat_tri_avg_var_4 mlat_tri_avg_neg_4 ];
x = categorical({'Front Delt', 'Longitudinal Tricep', 'Lateral Tricep'});
x = reordercats(x,{'Front Delt', 'Longitudinal Tricep', 'Lateral Tricep'});

bar(x,y); legend('Positive','Variable','Negative'); title("Average Muscle Activation for Movement Phase in Upwards Dir."); ylabel("Percent Activation");xlabel("Muscle");

% Downwards 
% Agonists are forearm, bicep, and rear-delt
subplot(1,2,2);
y = [ mfore_avg_pos_3 mfore_avg_var_3 mfore_avg_neg_3; mbi_avg_pos_3 mbi_avg_var_3 mbi_avg_neg_3; mrear_delt_avg_pos_3 mrear_delt_avg_var_3 mrear_delt_avg_neg_3];
x = categorical({'Forearm','Bicep','Rear Delt' });
x = reordercats(x,{'Forearm','Bicep','Rear Delt' });

bar(x,y); legend('Positive','Variable','Negative'); title("Average Muscle Activation for Movement Phase in Downwards Dir."); ylabel("Percent Activation");xlabel("Muscle");
%%

% Upwards 
figure; 
subplot (1,2,1);
y = [mfront_delt_avg_pos_4 mfront_delt_avg_var_4 mfront_delt_avg_neg_4; mlong_tri_avg_pos_4 mlong_tri_avg_var_4 mlong_tri_avg_neg_4; mlat_tri_avg_pos_4 mlat_tri_avg_var_4 mlat_tri_avg_neg_4 ];
x = categorical({'Front Delt', 'Longitudinal Tricep', 'Lateral Tricep'});
x = reordercats(x,{'Front Delt', 'Longitudinal Tricep', 'Lateral Tricep'});

bar(x,y); legend('Positive','Variable','Negative'); title("Average Muscle Activation for Movement Period in Upwards Dir."); ylabel("Percent Activation");xlabel("Muscle");

% Downwards 
% Agonists are forearm, bicep, and rear-delt
subplot(1,2,2);
y = [ mfore_avg_pos_3 mfore_avg_var_3 mfore_avg_neg_3; mbi_avg_pos_3 mbi_avg_var_3 mbi_avg_neg_3; mrear_delt_avg_pos_3 mrear_delt_avg_var_3 mrear_delt_avg_neg_3];
x = categorical({'Forearm','Bicep','Rear Delt' });
x = reordercats(x,{'Forearm','Bicep','Rear Delt' });

bar(x,y); legend('Positive','Variable','Negative'); title("Average Muscle Activation for Movement Period in Downwards Dir."); ylabel("Percent Activation");xlabel("Muscle");

%% Stability Phase 
% Upwards 
% Agonists are triceps and front delt
figure; 
subplot (1,2,1);
y = [mfront_delt_avg_pos_4 mfront_delt_avg_var_4 mfront_delt_avg_neg_4; mlong_tri_avg_pos_4 mlong_tri_avg_var_4 mlong_tri_avg_neg_4; mlat_tri_avg_pos_4 mlat_tri_avg_var_4 mlat_tri_avg_neg_4 ];
x = categorical({'Front Delt', 'Longitudinal Tricep', 'Lateral Tricep'});
x = reordercats(x,{'Front Delt', 'Longitudinal Tricep', 'Lateral Tricep'});

bar(x,y); legend('Positive','Variable','Negative'); title("Average Muscle Activation for Stability Phase in Upwards Dir."); ylabel("Percent Activation");xlabel("Muscle");

% Downwards 
% Agonists are forearm, bicep, and rear-delt
subplot(1,2,2);
y = [ mfore_avg_pos_3 mfore_avg_var_3 mfore_avg_neg_3; mbi_avg_pos_3 mbi_avg_var_3 mbi_avg_neg_3; mrear_delt_avg_pos_3 mrear_delt_avg_var_3 mrear_delt_avg_neg_3];
x = categorical({'Forearm','Bicep','Rear Delt' });
x = reordercats(x,{'Forearm','Bicep','Rear Delt' });

bar(x,y); legend('Positive','Variable','Negative'); title("Average Muscle Activation for Stability Phase in Downwards Dir."); ylabel("Percent Activation");xlabel("Muscle");

