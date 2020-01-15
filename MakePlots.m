%% Collect Trials from MATs

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
                'Subject28Data.mat'
              };

clear trials
trials = [];
for i = 1:length(subjectMats)
    file = fullfile(subjectMatFolder, subjectMats{i});
    temp = load(file, 'trials');
    trials = [trials; temp.trials];
end
          

% Rise Time

for i = 1:length(trials)
    
    x = trials(i).Data.EndEffPos_FromJA;
    
    zero_level = GetZeroLevel(x);
    ss_level = GetSteadyStateLevel(x);
    
    ten_percent_level = zero_level + 0.1*(ss_level - zero_level);
    ninty_percent_level = zero_level + 0.9*(ss_level - zero_level);
    
    if (ss_level < 0)
        inds = x < ten_percent_level;
        ten_percent_time = find(inds, 1, 'first');
        
        inds = x < ninty_percent_level;
        ninty_percent_time = find(inds, 1, 'first');
    else
        inds = x > ten_percent_level;
        ten_percent_time = find(inds, 1, 'first');
        
        inds = x > ninty_percent_level;
        ninty_percent_time = find(inds, 1, 'first');
        
    end
    trials(i).RiseTime = ninty_percent_time - ten_percent_time;
    
end

% Percent Overshoot
for i = 1:length(trials)
    x = trials(i).Data.EndEffPos_FromJA;

    trial_zero = mean(x(1:30));
    trial_ss = mean(x(end-100:end));
    trial_range = abs(trial_ss - trial_zero);
    trial_max = max(abs(x));

    trials(i).PercentOvershoot = (abs(trial_max - trial_zero)./trial_range - 1)*100;

end

%  Calculate settling time
for i = 1:length(trials)
    
    x = trials(i).Data.EndEffPos_FromJA;
    ss_val = GetSteadyStateLevel(x);
    zero_val = GetZeroLevel(x);
    x_range = abs(ss_val - zero_val);
        
    settleBound = 0.05*x_range;
    unsettledInds = (x < (ss_val - settleBound)) | ( x > (ss_val + settleBound));
    trials(i).SettlingTime = find(unsettledInds, 1, 'last');
end

% Energy
for i = 1:length(trials)
    x = trials(i).Data.EndEffPos_FromJA;
    fx = trials(i).Data.Force(:,1);
    dx = diff(x);
    d_energy = [0; fx(1:length(dx)).*dx];
    for iEnergy = 1:length(d_energy)
        energy(iEnergy) = sum(d_energy(1:iEnergy));
    end
    trials(i).Energy = energy;
    
    v = trials(i).Data.xdot;
    trials(i).Power = fx.*v;
    trials(i).PowerMax = max(trials(i).Power);
end



%% Make Plots Individually

figCount = 0;
for targetDirNum = 1:4
    % Make plots
    figs = [];
    for dampNum = 1:3
        figs = [figs, MakeX_VA_Damp_Plot(trials, targetDirNum, dampNum)];
    end

    % Fix Axes Limits
    
    % Get ylims for x and va subplots on all figures
    x_ylims = [];
    va_ylims = [];
    damp_ylims = [];
    for i = 1:length(figs)
        fig = figs(i);
        axes_all = fig.Children; % There should be 3 axes per fig b/c each fig has 3 subplots

        x_ylims = [x_ylims; ylim(axes_all(3))];
        va_ylims = [va_ylims; ylim(axes_all(2))];
        damp_ylims = [damp_ylims; ylim(axes_all(1))];
    end

    % Get the extreme y lims for x and va
    x_ylim = [min(x_ylims(:,1)), max(x_ylims(:,2))];
    va_ylim = [min(va_ylims(:,1)), max(va_ylims(:,2))];
    damp_ylim = [min(damp_ylims(:,1)), max(damp_ylims(:,2))];
    
    % Get range of limits 
    x_range = x_ylim(2) - x_ylim(1);
    va_range = va_ylim(2) - va_ylim(1);
    damp_range = damp_ylim(2) - damp_ylim(1);

    % Make limits 10% wider
    x_ylim(1) = x_ylim(1) - 0.1*x_range;
    x_ylim(2) = x_ylim(2) + 0.1*x_range;
    
    va_ylim(1) = va_ylim(1) - 0.1*va_range;
    va_ylim(2) = va_ylim(2) + 0.1*va_range;
    
    damp_ylim(1) = damp_ylim(1) - 0.1*damp_range;
    damp_ylim(2) = damp_ylim(2) + 0.1*damp_range;
    
    % Set ylims
    for i = 1:length(figs)
        fig = figs(i);
        axes_all = fig.Children;

        ylim(axes_all(1), damp_ylim);
        ylim(axes_all(2), va_ylim);
        ylim(axes_all(3), x_ylim);
    end
        
end

%% Make X, VA, Damp Plots Combined by Target Location (Left, Right, Down, Up)

locations = {'Left', 'Right', 'Down', 'Up'};
for targetDirNum = 1:4
    fig = Make_X_VA_Damp_Plot_Combined(trials, targetDirNum, 'ShowMovementMatch', false, 'ShowViolations', false);
%     set(fig,'PaperOrientation','landscape');
%     set(fig,'PaperUnits','normalized');
%     set(fig,'PaperUnits','inches','PaperPosition', [0 0 12 8]);
%     figName = sprintf('CombinedPlot_%s.jpg', locations{targetDirNum});
%     print(fig, '-djpeg', figName, '-r300')
    
end

%% Make max velocity bar graphs

subject = 22;

% Collect all filtered data
data_processed = [];
for targetNum = 1:4
    for dampNum = 1:3
        temp = FilterAndAnalyzeData(trials, targetNum, dampNum);
        data_processed = [data_processed; temp];
        
        % print data table
        maxVelAll = [];
        maxVelArr = temp.XdotMax;
        fprintf('\n\n')
        fprintf('Direction: %s        Damping: %s\n', temp.TargetDirText, temp.DampingText)    
        fprintf('--------------------------------------------------------------\n');
        fprintf('Subject\tMean[m/s]\tStd[m/s]\tMax[m/s]\tn\n')
        fprintf('--------------------------------------------------------------\n');
        for i =1:length(maxVelArr)
            maxVel = maxVelArr(i);
            fprintf('%d\t%f\t%f\t%f\t%d\n', ...
                    maxVel.SubjectNumber, ...
                    maxVel.max_xdot_mean, ...
                    maxVel.max_xdot_std, ...
                    maxVel.max_xdot_max, ...
                    maxVel.nTrials)
            maxVelAll = [maxVelAll, maxVel.max_xdot];
        end
        fprintf('--------------------------------------------------------------\n');
        fprintf('Total\t%f\t%f\t%f\t%d\n', ...
                mean(maxVelAll), ...
                std(maxVelAll), ...
                max(maxVelAll), ...
                length(maxVelAll));
        
        
        
    end
end

orange = [255, 165, 0]/256;
dampColors = {'b', 'g', orange};


fig = figure;
set(fig,'Color',[1,1,1]);
for targetNum = 1:4
    % Get subset of data of based on target number
    set_inds = [data_processed.TargetDirNum] == targetNum;
    data_processed_set = data_processed(set_inds);
    
    % Get means, std, damp label text
    xdotmax_mean = zeros(1,3);
    xdotmax_std = zeros(1,3);
    nTrials = zeros(1,3);
    for i = 1:3
        xdotmax_data = data_processed_set(i).XdotMax;
        subjectInds = [xdotmax_data.SubjectNumber];
        subjectInd = subjectInds == subject;
        
        xdotmax_subject = xdotmax_data(subjectInd);
        xdotmax_mean(i) = xdotmax_subject.max_xdot_mean;
        xdotmax_std(i) = xdotmax_subject.max_xdot_std;
        nTrials(i) = xdotmax_subject.nTrials;
    end
    
    dampTexts = {};
    dampTextsCat = {};
    for i = 1:length(data_processed_set)
        dText = sprintf('%s (n=%d)', data_processed_set(i).DampingText, ...
                                     nTrials(i));
        dampTexts{i} = dText;
        dampTextsCat{i} = data_processed_set(i).DampingText;
    end
    c = categorical(dampTextsCat);
    
    % Plot 
    ax = subplot(1, 4, targetNum);
    for i = 1:length(xdotmax_mean)
        bar(c(i),xdotmax_mean(i), 'FaceColor', dampColors{i});
        hold on
    end
    
    for i = 1:length(xdotmax_mean)
        er = errorbar(c(i),xdotmax_mean(i), xdotmax_std(i));
        er.Color = [0,0,0];
        er.LineStyle = 'none';
    end
    
    % Make plot pretty
    if (targetNum ~= 2)
        title(data_processed_set(1).TargetDirText);
    else
        title(sprintf('Subject %d\n\n%s', subject, data_processed_set(1).TargetDirText))
    end
    if (targetNum == 1)
        ylabel('$\dot{x}$ [m/s]', 'Interpreter', 'latex');
    end
    box('off')
    set(ax, 'ticklength', [0,0]);
    legend(dampTexts);
    legend('boxoff')
    hold off
end

% Set Y Limits equal
child = fig.Children;
axes_arr = [];
for i = length(child):-1:1
    if (~isa(child(i), 'matlab.graphics.axis.Axes'))
        child(i) = [];
    end
end
SetYLimsEqual(child, 'BottomPadPercent', 0)


%% Make percent overshoot bar graphs
subject = 22;

% po_rm_meas = zeros(10,4,4);


% Collect all filtered data
data_processed = [];
for targetNum = 1:4
    for dampNum = 1:3
        temp = FilterAndAnalyzeData(trials, targetNum, dampNum);
        data_processed = [data_processed; temp];
        
%         po_rm_means_temp = zeros(10,2);
        
        % print data table
        poAll = [];
        poArr = temp.PercentOvershoot;
        fprintf('\n\n')
        fprintf('Direction: %s        Damping: %s\n', temp.TargetDirText, temp.DampingText)    
        fprintf('--------------------------------------------------------------\n');
        fprintf('Percent Overshoot\n')
        fprintf('Subject\tMean\tStd\tMax\tn\n')
        fprintf('--------------------------------------------------------------\n');
        for i =1:length(poArr)
            po = poArr(i);
            fprintf('%d\t%f\t%f\t%f\t%d\n', ...
                    po.SubjectNumber, ...
                    po.po_mean, ...
                    po.po_std, ...
                    po.po_max, ...
                    po.nTrials)
            poAll = [poAll, po.po];
            
            % put data into temp repeated measures anova mean array
            po_rm_means_temp(i,1) = po.SubjectNumber;
            po_rm_means_temp(i,2) = po.po_mean;
        end
        fprintf('--------------------------------------------------------------\n');
        fprintf('Total\t%f\t%f\t%f\t%d\n', ...
                mean(poAll), ...
                std(poAll), ...
                max(poAll), ...
                length(poAll));
            
        % put temp data into repeated measures anova mean array 
        po_rm_means_temp = sortrows(po_rm_means_temp);
        po_rm_meas(:,1,targetNum) = po_rm_means_temp(:,1);
        po_rm_meas(:,dampNum+1,targetNum) = po_rm_means_temp(:,2);
            
    end
end

orange = [255, 165, 0]/256;
dampColors = {'b', 'g', orange};

fig = figure;
set(fig,'Color',[1,1,1]);
for targetNum = 1:4
    % Get subset of data of based on target number
    set_inds = [data_processed.TargetDirNum] == targetNum;
    data_processed_set = data_processed(set_inds);
    
    % Get means, std, damp label text
    po_mean = zeros(1,3);
    po_std = zeros(1,3);
    nTrials = zeros(1,3);
    for i = 1:3
        po_data = data_processed_set(i).PercentOvershoot;
        subjectInds = [po_data.SubjectNumber];
        subjectInd = subjectInds == subject;
        
        po_subject = po_data(subjectInd);
        po_mean(i) = po_subject.po_mean;
        po_std(i) = po_subject.po_std;
        nTrials(i) = po_subject.nTrials;
    end
    
    dampTexts = {};
    dampTextsCat = {};
    for i = 1:length(data_processed_set)
        dText = sprintf('%s (n=%d)', data_processed_set(i).DampingText, ...
                                     nTrials(i));
        dampTexts{i} = dText;
        dampTextsCat{i} = data_processed_set(i).DampingText;
    end
    c = categorical(dampTextsCat);
    
    % Plot 
    ax = subplot(1, 4, targetNum);
    for i = 1:length(po_mean)
        bar(c(i),po_mean(i), 'FaceColor', dampColors{i});
        hold on
    end
    
    for i = 1:length(po_mean)
        er = errorbar(c(i),po_mean(i), po_std(i));
        er.Color = [0,0,0];
        er.LineStyle = 'none';
    end
    
    % Make plot pretty
    if (targetNum ~= 2)
        title(data_processed_set(1).TargetDirText);
    else
        title(sprintf('Subject %d\n\n%s', subject, data_processed_set(1).TargetDirText))
    end
    if (targetNum == 1)
        ylabel('Percent Overshoot', 'Interpreter', 'latex');
    end
    box('off')
    set(ax, 'ticklength', [0,0]);
    legend(dampTexts);
    legend('boxoff')
    hold off
end

% Set Y Limits equal
child = fig.Children;
axes_arr = [];
for i = length(child):-1:1
    if (~isa(child(i), 'matlab.graphics.axis.Axes'))
        child(i) = [];
    end
end
SetYLimsEqual(child, 'BottomPadPercent', 0)

set(fig,'PaperUnits','inches','PaperPosition', [0 0 12 8]);
figName = sprintf('CombinedPlot_%s.jpg', locations{targetDirNum});
print(fig, '-djpeg', figName, '-r300')


%% Repeated Measures Anova - Percent Overshoot
fprintf('RM ANOVA - Percent Overshoot - All Damping Groups')
% Left
fprintf('\n\n\nLeft\n')
fprintf('Group1\tGroup2\tGroup3\tDF_{cond}\tDF_{error}\tF\tp\n');
res = CalcRepeatedAnova(po_rm_meas(:,2:4,1));
fprintf('%s\t%s\t%s\t%d\t%d\t%f\t%f\n', 'Pos.', 'Neg.', 'Var', res.DFcond, res.DFerror, res.F, res.p)

res = CalcRepeatedAnova(po_rm_meas(:,[2,3],1));
fprintf('%s\t%s\t%s\t%d\t%d\t%f\t%f\n', 'Pos.', 'Neg.', '', res.DFcond, res.DFerror, res.F, res.p)

res = CalcRepeatedAnova(po_rm_meas(:,[2,4],1));
fprintf('%s\t%s\t%s\t%d\t%d\t%f\t%f\n', 'Pos.', 'Var.', '', res.DFcond, res.DFerror, res.F, res.p)

res = CalcRepeatedAnova(po_rm_meas(:,[3,4],1));
fprintf('%s\t%s\t%s\t%d\t%d\t%f\t%f\n', 'Neg.', 'Var.', '', res.DFcond, res.DFerror, res.F, res.p)



% Right
fprintf('\n\n\nRight\n')
fprintf('Group1\tGroup2\tGroup3\tDF_{cond}\tDF_{error}\tF\tp\n');
res = CalcRepeatedAnova(po_rm_meas(:,2:4,2));
fprintf('%s\t%s\t%s\t%d\t%d\t%f\t%f\n', 'Pos.', 'Neg.', 'Var', res.DFcond, res.DFerror, res.F, res.p)

res = CalcRepeatedAnova(po_rm_meas(:,[2,3],2));
fprintf('%s\t%s\t%s\t%d\t%d\t%f\t%f\n', 'Pos.', 'Neg.', '', res.DFcond, res.DFerror, res.F, res.p)

res = CalcRepeatedAnova(po_rm_meas(:,[2,4],2));
fprintf('%s\t%s\t%s\t%d\t%d\t%f\t%f\n', 'Pos.', 'Var.', '', res.DFcond, res.DFerror, res.F, res.p)

res = CalcRepeatedAnova(po_rm_meas(:,[3,4],2));
fprintf('%s\t%s\t%s\t%d\t%d\t%f\t%f\n', 'Neg.', 'Var.', '', res.DFcond, res.DFerror, res.F, res.p)



% Down
fprintf('\n\n\nDown\n')
fprintf('Group1\tGroup2\tGroup3\tDF_{cond}\tDF_{error}\tF\tp\n');
res = CalcRepeatedAnova(po_rm_meas(:,2:4,3));
fprintf('%s\t%s\t%s\t%d\t%d\t%f\t%f\n', 'Pos.', 'Neg.', 'Var', res.DFcond, res.DFerror, res.F, res.p)

res = CalcRepeatedAnova(po_rm_meas(:,[2,3],3));
fprintf('%s\t%s\t%s\t%d\t%d\t%f\t%f\n', 'Pos.', 'Neg.', '', res.DFcond, res.DFerror, res.F, res.p)

res = CalcRepeatedAnova(po_rm_meas(:,[2,4],3));
fprintf('%s\t%s\t%s\t%d\t%d\t%f\t%f\n', 'Pos.', 'Var.', '', res.DFcond, res.DFerror, res.F, res.p)

res = CalcRepeatedAnova(po_rm_meas(:,[3,4],3));
fprintf('%s\t%s\t%s\t%d\t%d\t%f\t%f\n', 'Neg.', 'Var.', '', res.DFcond, res.DFerror, res.F, res.p)


% Up
fprintf('\n\n\nUp\n')
fprintf('Group1\tGroup2\tGroup3\tDF_{cond}\tDF_{error}\tF\tp\n');
res = CalcRepeatedAnova(po_rm_meas(:,2:4,4));
fprintf('%s\t%s\t%s\t%d\t%d\t%f\t%f\n', 'Pos.', 'Neg.', 'Var', res.DFcond, res.DFerror, res.F, res.p)

res = CalcRepeatedAnova(po_rm_meas(:,[2,3],4));
fprintf('%s\t%s\t%s\t%d\t%d\t%f\t%f\n', 'Pos.', 'Neg.', '', res.DFcond, res.DFerror, res.F, res.p)

res = CalcRepeatedAnova(po_rm_meas(:,[2,4],4));
fprintf('%s\t%s\t%s\t%d\t%d\t%f\t%f\n', 'Pos.', 'Var.', '', res.DFcond, res.DFerror, res.F, res.p)

res = CalcRepeatedAnova(po_rm_meas(:,[3,4],4));
fprintf('%s\t%s\t%s\t%d\t%d\t%f\t%f\n', 'Neg.', 'Var.', '', res.DFcond, res.DFerror, res.F, res.p)





%% Repeated Measures Anova - Percent Overshoot - Positive vs Negative
% Left
po_rm_results = [];
targetDirTexts = {'Left', 'Right', 'Down', 'Up'};
for targetNum = 1:4
    fprintf('\n\n\n%s\n', targetDirTexts{targetNum});
    fprintf('Positive vs Negative\n');
    temp = CalcRepeatedAnova(po_rm_meas(:,[2,3],targetNum));
    temp.Direction = targetDirTexts{targetNum};
    temp.DampGroup1 = 'Positive';
    temp.DampGroup2 = 'Negative';
    po_rm_results = [po_rm_results, temp];
    
    fprintf('\nPositive vs Variable\n');
    CalcRepeatedAnova(po_rm_meas(:,[2,4],targetNum));
    temp = CalcRepeatedAnova(po_rm_meas(:,[2,4],targetNum));
    temp.Direction = targetDirTexts{targetNum};
    temp.DampGroup1 = 'Positive';
    temp.DampGroup2 = 'Variable';
    po_rm_results = [po_rm_results, temp];

    fprintf('\nNegative vs Variable\n');
    temp = CalcRepeatedAnova(po_rm_meas(:,[3,4],targetNum));
    temp.Direction = targetDirTexts{targetNum};
    temp.DampGroup1 = 'Negative';
    temp.DampGroup2 = 'Variable';
    po_rm_results = [po_rm_results, temp];
end

%% Make Force Plots During Movement


% Collect all filtered data
data_processed = [];
for targetNum = 1:4
    for dampNum = 1:3
        data_processed = [data_processed; FilterAndAnalyzeData(trials, targetNum, dampNum)];
    end
end

% Line Color RGB vals
lightGreyColor = [179, 179, 179]/256;

fprintf('\n\n\n')
fprintf('--------------------------------------------------------------\n');
fprintf('Direction\tDamping\tSubject\tMeanMean[N]\tMeanStd[N]\tMaxMean[N]\tMaxStd[N]\tn\n')
fprintf('--------------------------------------------------------------\n');

x_all = {[],[],[];
         [],[],[];
         [],[],[];
         [],[],[]};
     
force_all = {[],[],[];
             [],[],[];
             [],[],[];
             [],[],[]};

force_rm_max = zeros(10,4,4);
force_rm_mean = zeros(10,4,4);
         
% Print force data
for targetNum = 1:4
    % Get subset of data of based on target number
    set_inds = [data_processed.TargetDirNum] == targetNum;
    data_processed_set = data_processed(set_inds);

    % Get all trials without violations with this target number
    trials_wov = [];
%     va_mean = {};
    va_mean_movement_length = [];
    for i = 1:length(data_processed_set)
        trials_wov = [trials_wov; data_processed_set(i).wov.trials];
        va_mean = data_processed_set(i).wov.va.mean_;
        va_mean_movement_length(end+1) = GetLastMoveInd(va_mean) - GetFirstMoveInd(va_mean);
    end
    
    max_movement_length = max(va_mean_movement_length);
    

    
    subjectNumbersAll = [trials_wov.SubjectNumber];
    subjectNumbers = unique(subjectNumbersAll);
    for i = 1:length(subjectNumbers)
        subjectNumber = subjectNumbers(i);
        subjectInds = subjectNumbersAll == subjectNumber;
        trials_subject = trials_wov(subjectInds);
        
        
 
        
        for dampNum = 1:3
            dampNumAll = [trials_subject.DampingNumber];
            dampNumInds = dampNumAll == dampNum;
            trials_subject_damp = trials_subject(dampNumInds);
            dampText = trials_subject_damp(i).DampingText;
            tarDirText = trials_subject_damp(i).TargetDirText;
            
            force = [];
            x = [];
            for j = 1:length(trials_subject_damp)
                firstMoveInd = trials_subject_damp(j).FirstMoveInd;
                lastMoveInd = firstMoveInd + max_movement_length;
                force = [force, trials_subject_damp(j).Data.Force(firstMoveInd:lastMoveInd, 1)];
                x = [x, trials_subject_damp(j).Data.EndEffPos_FromJA(firstMoveInd:lastMoveInd)];
            end
            
            % add raw data to arrays
            force_all{targetNum, dampNum} = [force_all{targetNum, dampNum}, force];
            x_all{targetNum, dampNum} = [x_all{targetNum, dampNum}, x];
            
            % calc mean, std, max of rms of the mean
            force_rms = sqrt(force.^2);
            force_rms_mean = mean(force_rms, 1);
            force_rms_max = max(force_rms, [], 1);
            
            force_rms_mean_mean = mean(force_rms_mean);
            force_rms_mean_std = std(force_rms_mean);
            force_rms_max_mean = mean(force_rms_max);
            force_rms_max_std = std(force_rms_max);
            
            % put data into repeated measures anova data array
            force_rm_mean(i,1,targetNum) = subjectNumber;
            force_rm_mean(i,dampNum+1,targetNum) = force_rms_mean_mean;
            force_rm_max(i,1,targetNum) = subjectNumber;
            force_rm_max(i,dampNum+1,targetNum) = force_rms_max_mean;
            
            % print data table
            fprintf('%s\t%s\t%d\t%f\t%f\t%f\t%f\t%d\n', ...
                    tarDirText, ...
                    dampText, ...
                    subjectNumber, ...
                    force_rms_mean_mean, ...
                    force_rms_mean_std, ...
                    force_rms_max_mean, ...
                    force_rms_max_std, ...
                    length(trials_subject_damp))


        end
    end
end

dampText = {'Positive', 'Negative', 'Variable'};
targetDirText = {'Left', 'Right', 'Down', 'Up'};

for targetNum = 1:4
    fig = figure;
    set(fig,'Color',[1,1,1]);
    for dampNum = 1:3
        x_temp = GetDataSetStats(x_all{targetNum, dampNum});
        force_temp = GetDataSetStats(force_all{targetNum, dampNum});
        
        % Plot x
        subplot(2,3,dampNum)
        plot(x_temp.mean_, 'Color', 'k');
        hold on
        plot(x_temp.mean_ + x_temp.std_, 'Color', 'k', 'LineStyle', '--');
        plot(x_temp.mean_ - x_temp.std_, 'Color', 'k', 'LineStyle', '--');
        hold off
        if (dampNum == 2)
            title(sprintf('%s\n\n%s', targetDirText{targetNum}, dampText{dampNum}));
        else
            title(sprintf('%s', dampText{dampNum}));
        end
        if (dampNum == 1)
            ylabel('x [m]')
        end
        
        
        % Plot force
        subplot(2,3,dampNum+3)
        force_rms = sqrt(force_temp.mean_.^2);
        plot(force_temp.mean_, 'k');
        hold on
        plot(force_temp.mean_ + force_temp.std_, 'k--');
        plot(force_temp.mean_ - force_temp.std_, 'k--', 'HandleVisibility', 'off');
        plot(mean(force_rms)*ones(size(force_rms)), 'r')
        plot(max(force_rms)*ones(size(force_rms)), 'b')
        hold on
        legend('Nominal', '+/-1 std', 'RMS Avg', '|Max|');
        legend('boxoff');
        legend('Location', 'southeast')
        if (dampNum == 1)
            ylabel('Force [N]')
        end
        xlabel('Time [ms]')
    end
    
    axes_arr = [subplot(2,3,1), subplot(2,3,2), subplot(2,3,3)];
    SetYLimsEqual(axes_arr);
    
    axes_arr = [subplot(2,3,4), subplot(2,3,5), subplot(2,3,6)];
    SetYLimsEqual(axes_arr);
    
    set(fig,'PaperUnits','inches','PaperPosition', [0 0 12 8]);
    figName = sprintf('Force_%s.jpg', targetDirText{targetNum});
    print(fig, '-djpeg', figName, '-r300')
    
end




fprintf('\n\n\n')
fprintf('--------------------------------------------------------------\n');
fprintf('Direction\tDamping\tMean[N]\tMax[N]\tn\n')
fprintf('--------------------------------------------------------------\n');

for targetNum = 1:4
    for dampNum = 1:3
        force_temp = GetDataSetStats(force_all{targetNum, dampNum});
        force_rms = sqrt(force_temp.mean_.^2);
        fprintf('%s\t%s\t%f\t%f\t%d\n', ...
                targetDirText{targetNum}, ...
                dampText{dampNum}, ...
                mean(force_rms), ...
                max(force_rms), ...
                size(force_temp.data_, 2));
                
        
    end
end

%%
% Repeated Measures Anova - Force Mean
% Left
fprintf('\n\n\nLeft - Force Mean\n')
CalcRepeatedAnova(force_rm_mean(:,2:4,1));

fprintf('\n\n\nLeft - Force Mean\n')
fprintf('Subject\tPositive\tNegative\tVariable\n')
for i = 1:size(force_rm_mean,1)
    fprintf('%d\t%f\t%f\t%f\n', force_rm_mean(i,1,1), ...
                                force_rm_mean(i,2,1), ...
                                force_rm_mean(i,3,1), ...
                                force_rm_mean(i,4,1))
end

% Right
fprintf('\n\n\nRight - Force Mean\n')
CalcRepeatedAnova(force_rm_mean(:,2:4,2));

fprintf('\n\n\nRight - Force Mean\n')
fprintf('Subject\tPositive\tNegative\tVariable\n')
for i = 1:size(force_rm_mean,1)
    fprintf('%d\t%f\t%f\t%f\n', force_rm_mean(i,1,2), ...
                                force_rm_mean(i,2,2), ...
                                force_rm_mean(i,3,2), ...
                                force_rm_mean(i,4,2))
end

% Down
fprintf('\n\n\nDown - Force Mean\n')
CalcRepeatedAnova(force_rm_mean(:,2:4,3));

fprintf('\n\n\nDown - Force Mean\n')
fprintf('Subject\tPositive\tNegative\tVariable\n')
for i = 1:size(force_rm_mean,1)
    fprintf('%d\t%f\t%f\t%f\n', force_rm_mean(i,1,3), ...
                                force_rm_mean(i,2,3), ...
                                force_rm_mean(i,3,3), ...
                                force_rm_mean(i,4,3))
end

% Up
fprintf('\n\n\nUp - Force Mean\n')
CalcRepeatedAnova(force_rm_mean(:,2:4,4));

fprintf('\n\n\nUp - Force Mean\n')
fprintf('Subject\tPositive\tNegative\tVariable\n')
for i = 1:size(force_rm_mean,1)
    fprintf('%d\t%f\t%f\t%f\n', force_rm_mean(i,1,4), ...
                                force_rm_mean(i,2,4), ...
                                force_rm_mean(i,3,4), ...
                                force_rm_mean(i,4,4))
end

% Repeated Measures Anova - Force Max
% Left
fprintf('\n\n\nLeft - Force Max\n')
CalcRepeatedAnova(force_rm_max(:,2:4,1));

fprintf('\n\n\nLeft - Force Max\n')
fprintf('Subject\tPositive\tNegative\tVariable\n')
for i = 1:size(force_rm_mean,1)
    fprintf('%d\t%f\t%f\t%f\n', force_rm_max(i,1,1), ...
                                force_rm_max(i,2,1), ...
                                force_rm_max(i,3,1), ...
                                force_rm_max(i,4,1))
end

% Right
fprintf('\n\n\nRight - Force Max\n')
CalcRepeatedAnova(force_rm_max(:,2:4,2));

fprintf('\n\n\nRight - Force Max\n')
fprintf('Subject\tPositive\tNegative\tVariable\n')
for i = 1:size(force_rm_mean,1)
    fprintf('%d\t%f\t%f\t%f\n', force_rm_max(i,1,2), ...
                                force_rm_max(i,2,2), ...
                                force_rm_max(i,3,2), ...
                                force_rm_max(i,4,2))
end

% Down
fprintf('\n\n\nDown - Force Max\n')
CalcRepeatedAnova(force_rm_max(:,2:4,3));

fprintf('\n\n\nDown - Force Max\n')
fprintf('Subject\tPositive\tNegative\tVariable\n')
for i = 1:size(force_rm_mean,1)
    fprintf('%d\t%f\t%f\t%f\n', force_rm_max(i,1,3), ...
                                force_rm_max(i,2,3), ...
                                force_rm_max(i,3,3), ...
                                force_rm_max(i,4,3))
end

% Up
fprintf('\n\n\nUp - Force Max\n')
CalcRepeatedAnova(force_rm_max(:,2:4,4));

fprintf('\n\n\nUp - Force Max\n')
fprintf('Subject\tPositive\tNegative\tVariable\n')
for i = 1:size(force_rm_mean,1)
    fprintf('%d\t%f\t%f\t%f\n', force_rm_max(i,1,4), ...
                                force_rm_max(i,2,4), ...
                                force_rm_max(i,3,4), ...
                                force_rm_max(i,4,4))
end


%% Force Max F p tables
% Left
fprintf('\n\n\nLeft\n')
fprintf('Group1\tGroup2\tGroup3\tDF_{cond}\tDF_{error}\tF\tp\n');
res = CalcRepeatedAnova(force_rm_max(:,2:4,1));
fprintf('%s\t%s\t%s\t%d\t%d\t%f\t%f\n', 'Pos.', 'Neg.', 'Var', res.DFcond, res.DFerror, res.F, res.p)

res = CalcRepeatedAnova(force_rm_max(:,[2,3],1));
fprintf('%s\t%s\t%s\t%d\t%d\t%f\t%f\n', 'Pos.', 'Neg.', '', res.DFcond, res.DFerror, res.F, res.p)

res = CalcRepeatedAnova(force_rm_max(:,[2,4],1));
fprintf('%s\t%s\t%s\t%d\t%d\t%f\t%f\n', 'Pos.', 'Var.', '', res.DFcond, res.DFerror, res.F, res.p)

res = CalcRepeatedAnova(force_rm_max(:,[3,4],1));
fprintf('%s\t%s\t%s\t%d\t%d\t%f\t%f\n', 'Neg.', 'Var.', '', res.DFcond, res.DFerror, res.F, res.p)



% Right
fprintf('\n\n\nRight\n')
fprintf('Group1\tGroup2\tGroup3\tDF_{cond}\tDF_{error}\tF\tp\n');
res = CalcRepeatedAnova(force_rm_max(:,2:4,2));
fprintf('%s\t%s\t%s\t%d\t%d\t%f\t%f\n', 'Pos.', 'Neg.', 'Var', res.DFcond, res.DFerror, res.F, res.p)

res = CalcRepeatedAnova(force_rm_max(:,[2,3],2));
fprintf('%s\t%s\t%s\t%d\t%d\t%f\t%f\n', 'Pos.', 'Neg.', '', res.DFcond, res.DFerror, res.F, res.p)

res = CalcRepeatedAnova(force_rm_max(:,[2,4],2));
fprintf('%s\t%s\t%s\t%d\t%d\t%f\t%f\n', 'Pos.', 'Var.', '', res.DFcond, res.DFerror, res.F, res.p)

res = CalcRepeatedAnova(force_rm_max(:,[3,4],2));
fprintf('%s\t%s\t%s\t%d\t%d\t%f\t%f\n', 'Neg.', 'Var.', '', res.DFcond, res.DFerror, res.F, res.p)



% Down
fprintf('\n\n\nDown\n')
fprintf('Group1\tGroup2\tGroup3\tDF_{cond}\tDF_{error}\tF\tp\n');
res = CalcRepeatedAnova(force_rm_max(:,2:4,3));
fprintf('%s\t%s\t%s\t%d\t%d\t%f\t%f\n', 'Pos.', 'Neg.', 'Var', res.DFcond, res.DFerror, res.F, res.p)

res = CalcRepeatedAnova(force_rm_max(:,[2,3],3));
fprintf('%s\t%s\t%s\t%d\t%d\t%f\t%f\n', 'Pos.', 'Neg.', '', res.DFcond, res.DFerror, res.F, res.p)

res = CalcRepeatedAnova(force_rm_max(:,[2,4],3));
fprintf('%s\t%s\t%s\t%d\t%d\t%f\t%f\n', 'Pos.', 'Var.', '', res.DFcond, res.DFerror, res.F, res.p)

res = CalcRepeatedAnova(force_rm_max(:,[3,4],3));
fprintf('%s\t%s\t%s\t%d\t%d\t%f\t%f\n', 'Neg.', 'Var.', '', res.DFcond, res.DFerror, res.F, res.p)


% Up
fprintf('\n\n\nUp\n')
fprintf('Group1\tGroup2\tGroup3\tDF_{cond}\tDF_{error}\tF\tp\n');
res = CalcRepeatedAnova(force_rm_max(:,2:4,4));
fprintf('%s\t%s\t%s\t%d\t%d\t%f\t%f\n', 'Pos.', 'Neg.', 'Var', res.DFcond, res.DFerror, res.F, res.p)

res = CalcRepeatedAnova(force_rm_max(:,[2,3],4));
fprintf('%s\t%s\t%s\t%d\t%d\t%f\t%f\n', 'Pos.', 'Neg.', '', res.DFcond, res.DFerror, res.F, res.p)

res = CalcRepeatedAnova(force_rm_max(:,[2,4],4));
fprintf('%s\t%s\t%s\t%d\t%d\t%f\t%f\n', 'Pos.', 'Var.', '', res.DFcond, res.DFerror, res.F, res.p)

res = CalcRepeatedAnova(force_rm_max(:,[3,4],4));
fprintf('%s\t%s\t%s\t%d\t%d\t%f\t%f\n', 'Neg.', 'Var.', '', res.DFcond, res.DFerror, res.F, res.p)



%% Rise Time


% Collect all trials without violations
trials_wov = [];
for targetNum = 1:4
    for dampNum = 1:3
        temp = FilterAndAnalyzeData(trials, targetNum, dampNum);
        trials_wov = [trials_wov; temp.wov.trials];
    end
end


rise_time_all = {[],[],[];
                [],[],[];
                [],[],[];
                [],[],[]};
dampText = {'Positive', 'Negative', 'Variable'};
targetDirText = {'Left', 'Right', 'Down', 'Up'};

rt_rm_meas = zeros(10, 4, 4);

fprintf('Rise Time\n')
fprintf('Direction\tDamping\tSubject\tMean[ms]\tStd[ms]\tn\n')

subjectNumbersAll = [trials_wov.SubjectNumber];
subjectNumbers = unique(subjectNumbersAll);
    
for iSubject = 1:length(subjectNumbers)
    subjectNumber = subjectNumbers(iSubject);

    subjectInds = subjectNumbersAll == subjectNumber;
    trials_subject = trials_wov(subjectInds);

    targetNumbers = [trials_subject.TargetDirNum];
    dampNumbers = [trials_subject.DampingNumber];

    rt_rm_meas(iSubject, 1, :) = subjectNumber;
    
    
    for targetNum = 1:4
        for dampNum = 1:3
            targetInds = targetNumbers == targetNum;
            dampInds = dampNumbers == dampNum;

            setInds = targetInds & dampInds;
            trials_set = trials_subject(setInds);

            rise_time_subject = [trials_set.RiseTime];
            
            if (dampNum == 3)
                mer = 1;
            end
            
            rt_rm_meas(iSubject, dampNum+1, targetNum) = mean(rise_time_subject);


            fprintf('%s\t%s\t%d\t%f\t%f\t%d\n', ...
                    targetDirText{targetNum}, ...
                    dampText{dampNum}, ...
                    subjectNumber, ...
                    mean(rise_time_subject), ...
                    std(rise_time_subject), ...
                    length(rise_time_subject));

            rise_time_all{targetNum, dampNum} = [rise_time_all{targetNum, dampNum}, rise_time_subject];

        end
    end

end

fprintf('\n\n')
fprintf('Rise Time\n')
fprintf('Direction\tDamping\tMean[ms]\tStd[ms]\tn\n')
for targetNum = 1:4
    for dampNum = 1:3

        rise_time = rise_time_all{targetNum, dampNum};

        fprintf('%s\t%s\t%f\t%f\t%d\n', ...
                targetDirText{targetNum}, ...
                dampText{dampNum}, ...
                mean(rise_time), ...
                std(rise_time), ...
                length(rise_time));

        rise_time_all{targetNum, dampNum} = [rise_time_all{targetNum, dampNum}, ];

    end
end

% Make plot
fig = figure;
set(fig,'Color',[1,1,1]);

c = categorical(dampText);
for targetNum = 1:4
    rise_time_means = zeros(1,3);
    rise_time_stds = zeros(1,3);
    for dampNum = 1:3
        rise_time_means(dampNum) = mean(rise_time_all{targetNum, dampNum});
        rise_time_stds(dampNum) = std(rise_time_all{targetNum, dampNum});
    end
    
    % Plot 
    ax = subplot(1, 4, targetNum);
    for i = 1:3
        bar(c(i),rise_time_means(i), 'FaceColor', dampColors{i});
        hold on
    end
    
    for i = 1:3
        er = errorbar(c(i),rise_time_means(i), rise_time_stds(i));
        er.Color = [0,0,0];
        er.LineStyle = 'none';
    end
    
    % Make plot pretty
    title(targetDirText{targetNum});
    if (targetNum == 1)
        ylabel('Rise Time[ms]', 'Interpreter', 'latex');
    end
    box('off')
    set(ax, 'ticklength', [0,0]);
%     legend(dampText);
%     legend('boxoff')
    hold off
    

end

SetYLimsEqual(fig.Children)

set(fig,'PaperUnits','inches','PaperPosition', [0 0 12 8]);
figName = sprintf('Rise_Time.jpg');
print(fig, '-djpeg', figName, '-r300')


% Rise Time Repeated Measures Anova

% Left
fprintf('\n\n\nLeft\n')
fprintf('Group1\tGroup2\tGroup3\tDF_{cond}\tDF_{error}\tF\tp\n');
res = CalcRepeatedAnova(rt_rm_meas(:,2:4,1));
fprintf('%s\t%s\t%s\t%d\t%d\t%f\t%f\n', 'Pos.', 'Neg.', 'Var', res.DFcond, res.DFerror, res.F, res.p)

res = CalcRepeatedAnova(rt_rm_meas(:,[2,3],1));
fprintf('%s\t%s\t%s\t%d\t%d\t%f\t%f\n', 'Pos.', 'Neg.', '', res.DFcond, res.DFerror, res.F, res.p)

res = CalcRepeatedAnova(rt_rm_meas(:,[2,4],1));
fprintf('%s\t%s\t%s\t%d\t%d\t%f\t%f\n', 'Pos.', 'Var.', '', res.DFcond, res.DFerror, res.F, res.p)

res = CalcRepeatedAnova(rt_rm_meas(:,[3,4],1));
fprintf('%s\t%s\t%s\t%d\t%d\t%f\t%f\n', 'Neg.', 'Var.', '', res.DFcond, res.DFerror, res.F, res.p)



% Right
fprintf('\n\n\nRight\n')
fprintf('Group1\tGroup2\tGroup3\tDF_{cond}\tDF_{error}\tF\tp\n');
res = CalcRepeatedAnova(rt_rm_meas(:,2:4,2));
fprintf('%s\t%s\t%s\t%d\t%d\t%f\t%f\n', 'Pos.', 'Neg.', 'Var', res.DFcond, res.DFerror, res.F, res.p)

res = CalcRepeatedAnova(rt_rm_meas(:,[2,3],2));
fprintf('%s\t%s\t%s\t%d\t%d\t%f\t%f\n', 'Pos.', 'Neg.', '', res.DFcond, res.DFerror, res.F, res.p)

res = CalcRepeatedAnova(rt_rm_meas(:,[2,4],2));
fprintf('%s\t%s\t%s\t%d\t%d\t%f\t%f\n', 'Pos.', 'Var.', '', res.DFcond, res.DFerror, res.F, res.p)

res = CalcRepeatedAnova(rt_rm_meas(:,[3,4],2));
fprintf('%s\t%s\t%s\t%d\t%d\t%f\t%f\n', 'Neg.', 'Var.', '', res.DFcond, res.DFerror, res.F, res.p)



% Down
fprintf('\n\n\nDown\n')
fprintf('Group1\tGroup2\tGroup3\tDF_{cond}\tDF_{error}\tF\tp\n');
res = CalcRepeatedAnova(rt_rm_meas(:,2:4,3));
fprintf('%s\t%s\t%s\t%d\t%d\t%f\t%f\n', 'Pos.', 'Neg.', 'Var', res.DFcond, res.DFerror, res.F, res.p)

res = CalcRepeatedAnova(rt_rm_meas(:,[2,3],3));
fprintf('%s\t%s\t%s\t%d\t%d\t%f\t%f\n', 'Pos.', 'Neg.', '', res.DFcond, res.DFerror, res.F, res.p)

res = CalcRepeatedAnova(rt_rm_meas(:,[2,4],3));
fprintf('%s\t%s\t%s\t%d\t%d\t%f\t%f\n', 'Pos.', 'Var.', '', res.DFcond, res.DFerror, res.F, res.p)

res = CalcRepeatedAnova(rt_rm_meas(:,[3,4],3));
fprintf('%s\t%s\t%s\t%d\t%d\t%f\t%f\n', 'Neg.', 'Var.', '', res.DFcond, res.DFerror, res.F, res.p)


% Up
fprintf('\n\n\nUp\n')
fprintf('Group1\tGroup2\tGroup3\tDF_{cond}\tDF_{error}\tF\tp\n');
res = CalcRepeatedAnova(rt_rm_meas(:,2:4,4));
fprintf('%s\t%s\t%s\t%d\t%d\t%f\t%f\n', 'Pos.', 'Neg.', 'Var', res.DFcond, res.DFerror, res.F, res.p)

res = CalcRepeatedAnova(rt_rm_meas(:,[2,3],4));
fprintf('%s\t%s\t%s\t%d\t%d\t%f\t%f\n', 'Pos.', 'Neg.', '', res.DFcond, res.DFerror, res.F, res.p)

res = CalcRepeatedAnova(rt_rm_meas(:,[2,4],4));
fprintf('%s\t%s\t%s\t%d\t%d\t%f\t%f\n', 'Pos.', 'Var.', '', res.DFcond, res.DFerror, res.F, res.p)

res = CalcRepeatedAnova(rt_rm_meas(:,[3,4],4));
fprintf('%s\t%s\t%s\t%d\t%d\t%f\t%f\n', 'Neg.', 'Var.', '', res.DFcond, res.DFerror, res.F, res.p)


%% Calculate Energy In a Given Trial
trialNum = 1000;

figure
ax1 = subplot(4,1,1);
x = trials(trialNum).Data.EndEffPos_FromJA;
fx = trials(trialNum).Data.Force(:,1);
plot(x)
ylabel('x [m]')
xlabel('Time [ms]')

dampText = trials(trialNum).DampingText;
dirText = trials(trialNum).TargetDirText;
title(sprintf('Direction: %s\tDamping: %s', dirText, dampText));

ax2 = subplot(4,1,2);
plot(fx)
ylabel('Force [N]')
xlabel('Time [ms]')

% Calc energy
dx = diff(x);
d_energy = [0; (fx(1:length(dx)).*dx)];

% add zero to beginning
energy_abs = zeros(size(d_energy));
energy = zeros(size(d_energy));
energy_only_pos = zeros(size(d_energy));

% only positive d_energy
d_energy_only_pos = (d_energy + abs(d_energy))./2;

for i = 1:length(d_energy)
    energy_abs(i) = sum(abs(d_energy(1:i)));
    energy(i) = sum(d_energy(1:i));
    energy_only_pos(i) = sum(d_energy_only_pos(1:i));
end


ax3 = subplot(4,1,3);

% plot(energy_abs)
% hold on 
plot(energy);
% plot(energy_only_pos);
% hold off
ylabel('Energy [Nm]')
xlabel('Time [ms]')
% legend(['\Sigma |fx * dx |'], ['\Sigma fx * dx'], ['\Sigma fx * dx of only'  char(10) 'positive energy contrib'])
% legend boxoff
% legend('Location', 'east')
% 
% figure


% Calculate power
dt = 1/1000;
power1 = d_energy/dt;

v = trials(trialNum).Data.xdot;
power2 = fx.*v;

subplot(4,1,4)
plot(power1)
hold on
plot(power2)
hold off
legend('p = dEnergy/dt', 'p = f*v');
xlabel('Time [ms]');
ylabel('Power [Nm/s]')

%% Print Percent Overshoot Data

trials_wov = [];
for targetDirNum = 1:4
    for dampNum = 1:3
        trials_wov = [trials_wov; GetTrialsWithoutViolations(trials, targetDirNum, dampNum)];
    end
end

targetDirNum_all = [trials_wov.TargetDirNum];
dampNum_all = [trials_wov.DampingNumber];
subjectNum_all = [trials_wov.SubjectNumber];
subjectNumList = sort(unique(subjectNum_all));

targetDirList = {'Left', 'Right', 'Down', 'Up'};
dampList = {'Positive', 'Negative', 'Variable'};

fprintf('Percent Overshoot\n');
fprintf('Direction\tDamping\tSubject\tMean\tStd\tn\n');
for targetDirNum = 1:4
    for dampNum = 1:3
        for iSubject = 1:length(subjectNumList)
            subjectNum = subjectNumList(iSubject);
            
            % Get Set of Trials correponding to target direction, damping,
            % and subject number
            targetNum_inds = (targetDirNum_all == targetDirNum);
            dampNum_inds = (dampNum_all == dampNum);
            subjectNum_inds = (subjectNum_all == subjectNum);
            
            set_inds = (targetNum_inds & dampNum_inds & subjectNum_inds);
            trials_set = trials_wov(set_inds);
            
            subject_rt = [trials_set.PercentOvershoot];
            fprintf('%s\t%s\t%d\t%f\t%f\t%d\n', ...
                    targetDirList{targetDirNum}, ...
                    dampList{dampNum}, ...
                    subjectNum, ...
                    mean(subject_rt), ...
                    std(subject_rt), ...
                    length(subject_rt))
        end
    end
end

%% Print Rise Time Data

trials_wov = [];
for targetDirNum = 1:4
    for dampNum = 1:3
        trials_wov = [trials_wov; GetTrialsWithoutViolations(trials, targetDirNum, dampNum)];
    end
end

targetDirNum_all = [trials_wov.TargetDirNum];
dampNum_all = [trials_wov.DampingNumber];
subjectNum_all = [trials_wov.SubjectNumber];
subjectNumList = sort(unique(subjectNum_all));

targetDirList = {'Left', 'Right', 'Down', 'Up'};
dampList = {'Positive', 'Negative', 'Variable'};

fprintf('Rise Time [ms]\n');
fprintf('Direction\tDamping\tSubject\tMean\tStd\tn\n');
for targetDirNum = 1:4
    for dampNum = 1:3
        for iSubject = 1:length(subjectNumList)
            subjectNum = subjectNumList(iSubject);
            
            % Get Set of Trials correponding to target direction, damping,
            % and subject number
            targetNum_inds = (targetDirNum_all == targetDirNum);
            dampNum_inds = (dampNum_all == dampNum);
            subjectNum_inds = (subjectNum_all == subjectNum);
            
            set_inds = (targetNum_inds & dampNum_inds & subjectNum_inds);
            trials_set = trials_wov(set_inds);
            
            subject_rt = [trials_set.RiseTime];
            fprintf('%s\t%s\t%d\t%f\t%f\t%d\n', ...
                    targetDirList{targetDirNum}, ...
                    dampList{dampNum}, ...
                    subjectNum, ...
                    mean(subject_rt), ...
                    std(subject_rt), ...
                    length(subject_rt))
        end
    end
end


%% Print Settling Time Data

trials_wov = [];
for targetDirNum = 1:4
    for dampNum = 1:3
        trials_wov = [trials_wov; GetTrialsWithoutViolations(trials, targetDirNum, dampNum)];
    end
end

targetDirNum_all = [trials_wov.TargetDirNum];
dampNum_all = [trials_wov.DampingNumber];
subjectNum_all = [trials_wov.SubjectNumber];
subjectNumList = sort(unique(subjectNum_all));

targetDirList = {'Left', 'Right', 'Down', 'Up'};
dampList = {'Positive', 'Negative', 'Variable'};

fprintf('Settling Time [ms]\n');
fprintf('Direction\tDamping\tSubject\tMean\tStd\tn\n');
for targetDirNum = 1:4
    for dampNum = 1:3
        for iSubject = 1:length(subjectNumList)
            subjectNum = subjectNumList(iSubject);
            
            % Get Set of Trials correponding to target direction, damping,
            % and subject number
            targetNum_inds = (targetDirNum_all == targetDirNum);
            dampNum_inds = (dampNum_all == dampNum);
            subjectNum_inds = (subjectNum_all == subjectNum);
            
            set_inds = (targetNum_inds & dampNum_inds & subjectNum_inds);
            trials_set = trials_wov(set_inds);
            
            subject_st = [trials_set.SettlingTime];
            fprintf('%s\t%s\t%d\t%f\t%f\t%d\n', ...
                    targetDirList{targetDirNum}, ...
                    dampList{dampNum}, ...
                    subjectNum, ...
                    mean(subject_st), ...
                    std(subject_st), ...
                    length(subject_st))
            
        end
    end
end


for targetDirNum = 1:4
    figure
    st_means = zeros(1,3);
    for dampNum = 1:3
        % Get Set of Trials correponding to target direction, damping,
        % and subject number
        targetNum_inds = (targetDirNum_all == targetDirNum);
        dampNum_inds = (dampNum_all == dampNum);
        
        set_inds = (targetNum_inds & dampNum_inds);
        trials_set = trials_wov(set_inds);
        
        st_set = [trials_set.SettlingTime];
        st_mean(dampNum) = mean(st_set);
        dampTextsCat{dampNum} = trials_set(1).DampingText;
       
        % Histogram of settling time
        ax = subplot(1,3,dampNum)
        hist(st_set, 25);
        titleStr = sprintf('%s - %s', targetDirList{targetDirNum}, ...
                                      dampTextsCat{dampNum});
        title(titleStr);
        xlabel('Settling Time [ms]');
        ax.XLim = [0, 4000];
        
    end
    
    % average settling time
%     c = categorical(dampTextsCat);
%     
%     % Plot 
%     bar(c, st_mean);
    
end


%% Print Settling Time Data

trials_wov = [];
for targetDirNum = 1:4
    for dampNum = 1:3
        trials_wov = [trials_wov; GetTrialsWithoutViolations(trials, targetDirNum, dampNum)];
    end
end

targetDirNum_all = [trials_wov.TargetDirNum];
dampNum_all = [trials_wov.DampingNumber];
subjectNum_all = [trials_wov.SubjectNumber];
subjectNumList = sort(unique(subjectNum_all));

targetDirList = {'Left', 'Right', 'Down', 'Up'};
dampList = {'Positive', 'Negative', 'Variable'};

fprintf('Energy Inputted [Nm]\n');
fprintf('Direction\tDamping\tSubject\tMean\tStd\tn\n');
for targetDirNum = 1:4
    
    % find minimum length
    targetNum_inds = (targetDirNum_all == targetDirNum);
    trials_wov_dir = trials_wov(targetNum_inds);
    trial_lengths = [trials_wov_dir.nSamples];
    trial_length_dir = min(trial_lengths);
    
    for dampNum = 1:3
        for iSubject = 1:length(subjectNumList)
            subjectNum = subjectNumList(iSubject);
            
            % Get Set of Trials correponding to target direction, damping,
            % and subject number
            targetNum_inds = (targetDirNum_all == targetDirNum);
            dampNum_inds = (dampNum_all == dampNum);
            subjectNum_inds = (subjectNum_all == subjectNum);
            
            set_inds = (targetNum_inds & dampNum_inds & subjectNum_inds);
            trials_set = trials_wov(set_inds);
            
            % Get total energy inputted at the trial_length_dir sample
            subject_energy = zeros(1, length(trials_set));
            for iSubjectTrial = 1:length(trials_set)
                subject_energy(iSubjectTrial) = trials_set(iSubjectTrial).Energy(trial_length_dir);
            end
            fprintf('%s\t%s\t%d\t%f\t%f\t%d\n', ...
                    targetDirList{targetDirNum}, ...
                    dampList{dampNum}, ...
                    subjectNum, ...
                    mean(subject_energy), ...
                    std(subject_energy), ...
                    length(subject_energy))
            
        end
    end
end


% for targetDirNum = 1:4
%     figure
%     st_means = zeros(1,3);
%     for dampNum = 1:3
%         % Get Set of Trials correponding to target direction, damping,
%         % and subject number
%         targetNum_inds = (targetDirNum_all == targetDirNum);
%         dampNum_inds = (dampNum_all == dampNum);
%         
%         set_inds = (targetNum_inds & dampNum_inds);
%         trials_set = trials_wov(set_inds);
%         
%         st_set = [trials_set.SettlingTime];
%         st_mean(dampNum) = mean(st_set);
%         dampTextsCat{dampNum} = trials_set(1).DampingText;
%        
%         % Histogram of settling time
%         ax = subplot(1,3,dampNum)
%         hist(st_set, 25);
%         titleStr = sprintf('%s - %s', targetDirList{targetDirNum}, ...
%                                       dampTextsCat{dampNum});
%         title(titleStr);
%         ax.XLim = [0, 4000];
%         
%     end
%     
%     % average settling time
% %     c = categorical(dampTextsCat);
% %     
% %     % Plot 
% %     bar(c, st_mean);
%     
% end


%% Power (Max) - Print stats

trials_wov = [];
for targetDirNum = 1:4
    for dampNum = 1:3
        trials_wov = [trials_wov; GetTrialsWithoutViolations(trials, targetDirNum, dampNum)];
    end
end

trial_length = min([trials_wov.nSamples]);

targetDirNum_all = [trials_wov.TargetDirNum];
dampNum_all = [trials_wov.DampingNumber];
subjectNum_all = [trials_wov.SubjectNumber];
subjectNumList = sort(unique(subjectNum_all));

targetDirList = {'Left', 'Right', 'Down', 'Up'};
dampList = {'Positive', 'Negative', 'Variable'};

fprintf('Power Inputted [Nm/s]\n');
fprintf('Direction\tDamping\tSubject\tMean\tStd\tn\n');

subject_power_max_mean = zeros(4,3,10);
subject_power_max_stds = zeros(4,3,10);

subject_x_all = {};
subject_v_all = {};
subject_f_all = {};
subject_p_all = {};

for targetDirNum = 1:4
    
    % find minimum length
    targetNum_inds = (targetDirNum_all == targetDirNum);
    trials_wov_dir = trials_wov(targetNum_inds);
    trial_lengths = [trials_wov_dir.nSamples];
    trial_length_dir = min(trial_lengths);
    
    for dampNum = 1:3
        for iSubject = 1:length(subjectNumList)
            subjectNum = subjectNumList(iSubject);
            
            % Get Set of Trials correponding to target direction, damping,
            % and subject number
            targetNum_inds = (targetDirNum_all == targetDirNum);
            dampNum_inds = (dampNum_all == dampNum);
            subjectNum_inds = (subjectNum_all == subjectNum);
            
            set_inds = (targetNum_inds & dampNum_inds & subjectNum_inds);
            trials_set = trials_wov(set_inds);
            
            if (iSubject == 1 && targetDirNum == 3)
                nTrials = length(trials_set);
                subject_x = zeros(nTrials, trial_length);
                subject_v = zeros(nTrials, trial_length);
                subject_f = zeros(nTrials, trial_length);
                subject_p = zeros(nTrials, trial_length);
                for iTrial = 1:length(trials_set)
                    subject_x(iTrial,:) = trials_set(iTrial).Data.EndEffPos_FromJA(1:trial_length);
                    subject_v(iTrial,:) = trials_set(iTrial).Data.xdot(1:trial_length);
                    subject_f(iTrial,:) = trials_set(iTrial).Data.Force(1:trial_length);
                    subject_p(iTrial,:) = trials_set(iTrial).Power(1:trial_length);
                end
                
                subject_x_all{dampNum} = subject_x;
                subject_v_all{dampNum} = subject_v;
                subject_f_all{dampNum} = subject_f;
                subject_p_all{dampNum} = subject_p;
                
            end
            
            % Get total energy inputted at the trial_length_dir sample
            subject_powermax = [trials_set.PowerMax];
            fprintf('%s\t%s\t%d\t%f\t%f\t%d\n', ...
                    targetDirList{targetDirNum}, ...
                    dampList{dampNum}, ...
                    subjectNum, ...
                    mean(subject_powermax), ...
                    std(subject_powermax), ...
                    length(subject_powermax))
            
            subject_power_max_mean(targetDirNum, dampNum, iSubject) = mean(subject_powermax);
            subject_power_max_stds(targetDirNum, dampNum, iSubject) = std(subject_powermax);
                
        end
    end
end


orange = [255, 165, 0]/256;
dampColors = {'b', 'g', orange};

dampTexts = {'Positive', 'Negative', 'Variable'};
c = categorical(dampTexts);
subjectNum = 9;

% subject_power_max_mean = mean(subject_power_max_mean, 3);
% subject_power_max_stds = std(subject_power_max_stds,[],3);

fig = figure;
for targetDirNum = 1:4
    % Plot 
    ax = subplot(1, 4, targetDirNum);
    for dampNum = 1:3
        bar(c(dampNum),subject_power_max_mean(targetDirNum, dampNum, subjectNum), 'FaceColor', dampColors{dampNum});
        hold on

        er = errorbar(c(dampNum), subject_power_max_mean(targetDirNum, dampNum, subjectNum), ...
                            subject_power_max_stds(targetDirNum, dampNum, subjectNum));
        er.Color = [0,0,0];
        er.LineStyle = 'none';
        title(targetDirList{targetDirNum})
    end
        
end

subplot(1,4,1);
ylabel('Average Max Power [Nm/s] - Subject 29')



% Set Y Limits equal
child = fig.Children;
axes_arr = [];
for i = length(child):-1:1
    if (~isa(child(i), 'matlab.graphics.axis.Axes'))
        child(i) = [];
    end
end
SetYLimsEqual(child, 'BottomPadPercent', 0)


%% X - V - F - P plot
figure
for dampNum = 1:3
    subplot(4,1,1)
    plot(mean(subject_x_all{dampNum}), 'Color', dampColors{dampNum});
    ylabel('Position [m]')
    hold on

    subplot(4,1,2)
    plot(mean(subject_v_all{dampNum}), 'Color', dampColors{dampNum});
    ylabel('Velocity [m/s]')
    hold on

    subplot(4,1,3)
    plot(mean(subject_f_all{dampNum}), 'Color', dampColors{dampNum});
    ylabel('Force [N]')
    hold on

    subplot(4,1,4)
    plot(mean(subject_p_all{dampNum}), 'Color', dampColors{dampNum})
    ylabel('Power [Nm/s]')
    hold on
end
subplot(4,1,1)
title('Subject 21 - Down')

subplot(4,1,1)
xlabel('Time [ms]')
legend('Positive', 'Negative', 'Variable')
legend boxoff

for i = 1:4
    subplot(4,1,i)
    xlim([0,trial_length])
end