% Step 1: Run 1 of next 4 sections, which loads data into data_all (data copied
% straight from excel sheet).  

% Step 2: Run 'Make Plots' section, this will split data_all into means and
% standard error matrices, convert SE to CI, and make plot

%% Rise Time
title_part = 'Rise Time';
plot_name = 'rise_time';

data_all = [277.6	17.4	281.2	16.2	364.6	26.3	383.8	24.9	326.8	12.6;
            343	23.9	358.4	25.1	485.2	29.4	463.9	22.9	412.6	15.3;
            281.4	20.2	281.8	17.1	365.8	25.8	351	23	320	11.9];
        

%% Percent Overshoot
title_part = 'Percent Overshoot';
plot_name = 'percent_overshoot';

data_all = [11.2	2.1	13.1	2.2	11.3	1.7	16	2.4	12.9	1.1;
            3.7	0.6	4.1	0.5	4.4	0.4	4.8	0.7	4.2	0.3;
            6.3	0.9	7.1	1	5.3	0.5	7.2	1	6.5	0.4];

%% Max RMS Force
title_part = 'Max Interaction Force';
plot_name = 'max_force';

data_all = [28.3	2.8	25.8	2.1	16.9	1.5	17.1	2	22	1.3;
            34.5	3.1	32.1	3.1	21.7	2.2	22.1	2.3	27.6	1.5;
            23.1	2.2	21.3	1.9	13.9	1.7	14.3	1.9	18.2	1.1];


%% Mean RMS Force
title_part = 'Mean Interaction Force';
plot_name = 'mean_force';

data_all = [8.3	0.6	9.1	0.7	6.9	0.6	4.4	0.3	7.2	0.4;
            7.5	0.4	7.8	0.4	6.7	0.3	3.8	0.2	6.4	0.3;
            6.2	0.4	6.6	0.5	4.7	0.4	2.8	0.2	5.1	0.3];



%% Make Plots

targetDirText = {'Forward', 'Backward', 'Left', 'Right', 'All'};
[r, c] = size(data_all);
data_mean = [];
data_se = [];
for i = 1:c
    if (mod(i,2) == 1)
        data_mean(:,end+1) = data_all(:,i);
    else
        data_se(:,end+1) = data_all(:,i);
    end   
end

% Calc confidence interval from standard error
data_ci = data_se.*(1.96);

orange = [255, 165, 0]/256;
dampColors = {'g','b', orange};

dampTexts = {'Neg', 'Pos', 'Var'};
c = categorical(dampTexts);
subjectNum = 9;


fig = figure;
set(gcf,'Color',[1,1,1]);
for targetDirNum = 1:5
    % Plot 
    ax = subplot(1,5, targetDirNum);
    for dampNum = 1:3
        % Bar (Mean)
        bar(c(dampNum),data_mean(dampNum, targetDirNum), 'FaceColor', dampColors{dampNum}, 'EdgeColor', 'none');
        hold on

        % Errorbar (Confidence Interval)
        er = errorbar(c(dampNum), data_mean(dampNum, targetDirNum), ...
                                  data_ci(dampNum, targetDirNum));
        er.Color = [0,0,0];
        er.LineStyle = 'none';
        
        % Title
        if (targetDirNum == 3)
            titlestr = sprintf('%s\n\n%s', title_part, targetDirText{targetDirNum});
        else
            titlestr = targetDirText{targetDirNum};
        end
        
        title(titlestr)
        
        % Y Label
        if (targetDirNum == 1)
            if strcmp(title_part, 'Rise Time')
                ylabelstr = 'Time [ms]';
            elseif (strcmp(title_part, 'Percent Overshoot'))
                ylabelstr = '% Overshoot';
            elseif (strcmp(title_part, 'Max Interaction Force'))
                ylabelstr = 'Force [N]';
            elseif (strcmp(title_part, 'Mean Interaction Force'))
                ylabelstr = 'Force [N]';
            end
            ylabel(ylabelstr)
        end
    end
    % Turn axes box off
    set(ax, 'box','off')
end

% Set Y Limits equal
child = fig.Children;
axes_arr = [];
for i = length(child):-1:1
    if (~isa(child(i), 'matlab.graphics.axis.Axes'))
        child(i) = [];
    end
end

SetYLimsEqual(child, 'BottomPadPercent', 0, 'TopPadPercent', 0)

% Make RAL_Plots directory if necessary
if (exist('RAL_Plots', 'dir') ~= 7)
    mkdir('RAL_Plots');
end

% Save figure as png and fig files
png_file = strcat(plot_name, '.png');
fig_file = strcat(plot_name, '.fig');
png_file = fullfile(pwd, 'RAL_Plots', png_file);
fig_file = fullfile(pwd, 'RAL_Plots', fig_file);
saveas(fig, png_file);
saveas(fig, fig_file);

