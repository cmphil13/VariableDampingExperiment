%% Data Import

% Get files in subject directory
base_dir = uigetdir(pwd, 'Choose Subject Directory');
files_all = GetAllFilesInDir(base_dir);

% Get Kuka Data files
trials = [];
for i = 1:length(files_all)
    if (strcmp(files_all(i).name, 'KukaData.txt'))
        trials = [trials; files_all(i)];
    end
end
% 
% Calculate end effector x position as a function of joint angles
% Make KUKA robot with Corke Robotic Toolbox 
alpha = [pi/2, -pi/2, -pi/2, pi/2, pi/2, -pi/2, 0];
a = [0, 0, 0, 0, 0, 0, 0];
d = [0.36, 0, 0.42, 0, 0.4, 0, 0.126];
theta=[0 0 0 0 0 0 -0.958709];
dh = [theta' d' a' alpha'];

% BUILD ROBOT--------------------------------------------------------------
%\
for i = 1:length(dh(:,1))
    L{i} = Link('d', dh(i,2), 'a', dh(i,3), 'alpha', dh(i,4));
end
robot = SerialLink([L{1} L{2} L{3} L{4} L{5} L{6} L{7}]);


% Import data, calculate position
for i = length(trials):-1:1
    fprintf('Importing %d of %d\n', (length(trials) - i +1), length(trials));
    
    % Split folder path parts into cell array
    folder = trials(i).folder;
    if (ismac || isunix)
        filesep = '/';
    elseif (ispc)
        filesep = '\'; 
    else
        error('Platform not supported\n');
    end
	folder_parts = regexp(folder, filesep, 'split');

    % Get subject number, etc.
    % Assumed directory structure:
    %   /somepath/Suject{#}/{MovementDirection}/Group{#}/Trial{#}
    %
    % This will translate to:
    %   folder_parts{1:end-4}   = parts of somepath
    %   folder_part{end-3}      = 'Subject{#}'
    %   folder_part{end-2}      = '{MovementDirection}' --> This will be either 'Down_up' or 'Left_Right'
    %   folder_part{end-1}      = 'Group{#}'
    %   folder_part{end}        = 'Trial{#}'
    % Subject Number
    subjectStr = 'Subject';
    if (strcmp(folder_parts{end-3}(1:length(subjectStr)), subjectStr))
        trials(i).SubjectNumber = str2num(folder_parts{end-3}(length(subjectStr)+1:end));
    else
        fprintf('File: %s\n', fullfile(trials(i).folder, trials(i).name));
        fprintf('Subject Number not detected.  File removed from file list\n');
        trials(i) = [];
        continue;
    end
    
    % Movement direction
    lrStr = 'Left_Right';
    duStr = 'Down_Up';
    if (strcmp(folder_parts{end-2}, lrStr))
        trials(i).MoveDir = 1;
        trials(i).MoveDirText = 'Left-Right';
    elseif (strcmp(folder_parts{end-2}, duStr))
        trials(i).MoveDir = 2;
        trials(i).MoveDirText = 'Down-Up';
    else
        fprintf('File: %s\n', fullfile(trials(i).folder, trials(i).name));
        fprintf('Movement Direction not detected.  File removed from file list\n');
        trials(i) = [];
        continue;
    end
    
    % Group Number
    groupStr = 'Group';
    if (strcmp(folder_parts{end-1}(1:length(groupStr)), groupStr))
        trials(i).GroupNumber = str2num(folder_parts{end-1}(length(groupStr)+1:end));
    else
        fprintf('File: %s\n', fullfile(trials(i).folder, trials(i).name));
        fprintf('Group Number not detected.  File removed from file list\n');
        trials(i) = [];
        continue;
    end
    
    % Trial Number
    trialStr = 'Trial';
    if (strcmp(folder_parts{end}(1:length(trialStr)), trialStr))
        trials(i).TrialNumber = str2num(folder_parts{end}(length(trialStr)+1:end));
    else
        fprintf('File: %s\n', fullfile(trials(i).folder, trials(i).name));
        fprintf('Trial Number not detected.  File removed from file list\n');
        trials(i) = [];
        continue;
    end 
    
    % Load Kuka Data
    trials(i).Data = LoadKukaDataFile(fullfile(trials(i).folder, trials(i).name));

    % Load EMG Data (if it exists)
    trials(i).EmgFileName = 'EmgData.txt';
    trials(i).Data.Emg = LoadEmgDataFile(fullfile(trials(i).folder, trials(i).EmgFileName));
    if (numel(trials(i).Data.Emg) > 0)
        trials(i).HasEmgData = true;
    else
        trials(i).HasEmgData = false;
    end
    
    
    % Damping
    damp = trials(i).Data.Damping;
    if (std(damp) < 1e-6)   % Constant Damping
        if (damp(1) > 0)     % Positive Damping
            trials(i).DampingNumber = 1;
            trials(i).DampingText = 'Positive';
        elseif (damp(1) < 0) % Negative Damping
            trials(i).DampingNumber = 2;
            trials(i).DampingText = 'Negative';
        elseif (damp(1) == 0) % Zero Damping
            trials(i).DampingNumber = 4;
            trials(i).DampingText = 'Zero';
        end
    else                    % Variable Damping
        trials(i).DampingNumber = 3;
        trials(i).DampingText = 'Variable';
    end
    
    % Target Direction
    trials(i).TargetDirNum = trials(i).Data.TargetDirection(1);
    if (trials(i).TargetDirNum == 1)      % Left
        trials(i).TargetDirText = 'Left';
    elseif (trials(i).TargetDirNum == 2)    % Right
        trials(i).TargetDirText = 'Right';
    elseif (trials(i).TargetDirNum == 3)    % Down
        trials(i).TargetDirText = 'Down';
    elseif (trials(i).TargetDirNum == 4)    % Up
        trials(i).TargetDirText = 'Up';
    end

    % Combine damping type and trial direction into single number
    nTargetDirections = 4;
    directionNum = trials(i).TargetDirNum;
    dampNum = trials(i).DampingNumber;
    trials(i).DampDirNum = nTargetDirections*(directionNum - 1) + dampNum;
    
    % Get number of samples in trial
    trials(i).nSamples = size(trials(i).Data.JointAngle, 1);   

    % Calculate end effector position based on measured joint angles
    ja = trials(i).Data.JointAngle; % get joint angles
    T = robot.fkine(ja);
    x = T(1,4,:);                       % x coordinates of end-effector
    x = reshape(x,size(x,3), 1);
    trials(i).Data.EndEffPos_FromJA = x;
    
    % Calcualte velocity-acceleration parameter
    v = trials(i).Data.xdot;
    a = trials(i).Data.xdotdot;
    trials(i).Data.va = v.*a;
    
    % Calculate first and last movement indices
    trials(i).FirstMoveInd = GetFirstMoveInd(trials(i).Data.va);
    trials(i).LastMoveInd = GetLastMoveInd(trials(i).Data.va);
    
end

%% Combine to Struct
matfilename = sprintf('Subject%dData.mat', trials(1).SubjectNumber);
save(matfilename, 'trials');


%% Get Kp and Kn from Groups 0 and 1 - Left/Right Trials

% Get subset of trials with left/right movement
moveDirNums = [trials.MoveDir];
md_ind = moveDirNums == 1;
trials_md = trials(md_ind);

% Get subset of trials with group number 0 or 1
groupNumbers = [trials_md.GroupNumber];
gn0 = groupNumbers == 0;
gn1 = groupNumbers == 1;
gn01 = gn0 | gn1;

trials_gn01 = trials_md(gn01);
 
% Calculate average Kp, Kn for trials
Kp_arr = [];
Kn_arr = [];
for i = 1:length(trials_gn01)
    [Kp, Kn] = GetKs(trials_gn01(i).Data.xdot, trials_gn01(i).Data.xdotdot);
    Kp_arr = [Kp_arr, Kp];
    Kn_arr = [Kn_arr, Kn];
end

fprintf("Groups 0 & 1\n");
fprintf("Kp: %f\n", mean(Kp_arr));
fprintf("Kn: %f\n", mean(Kn_arr));

%% Get Kp and Kn from Groups 2 and 3  - Left/Right Trials

% Get subset of trials with left/right movement
moveDirNums = [trials.MoveDir];
md_ind = moveDirNums == 1;
trials_md = trials(md_ind);

% Get subset of trials with group number 0 or 1
groupNumbers = [trials_md.GroupNumber];
gn2 = groupNumbers == 2;
gn3 = groupNumbers == 3;
gn23 = gn2 | gn3;

trials_gn23 = trials_md(gn23);

% Calculate average Kp, Kn for trials
Kp_arr = [];
Kn_arr = [];
for i = 1:length(trials_gn23)
    [Kp, Kn] = GetKs(trials_gn23(i).Data.xdot, trials_gn23(i).Data.xdotdot);
    Kp_arr = [Kp_arr, Kp];
    Kn_arr = [Kn_arr, Kn];
end

fprintf("Groups 2 & 3\n");
fprintf("Kp: %f\n", mean(Kp_arr));
fprintf("Kn: %f\n", mean(Kn_arr));


%% Get Kp and Kn from Groups 0 and 1 - Down-Up Trials

% Get subset of trials with left/right movement
moveDirNums = [trials.MoveDir];
md_ind = moveDirNums == 2;
trials_md = trials(md_ind);

% Get subset of trials with group number 0 or 1
groupNumbers = [trials_md.GroupNumber];
gn0 = groupNumbers == 0;
gn1 = groupNumbers == 1;
gn01 = gn0 | gn1;

trials_gn01 = trials_md(gn01);

% Calculate average Kp, Kn for trials
Kp_arr = [];
Kn_arr = [];
for i = 1:length(trials_gn01)
    [Kp, Kn] = GetKs(trials_gn01(i).Data.xdot, trials_gn01(i).Data.xdotdot);
    Kp_arr = [Kp_arr, Kp];
    Kn_arr = [Kn_arr, Kn];
end

fprintf("Groups 0 & 1\n");
fprintf("Kp: %f\n", mean(Kp_arr));
fprintf("Kn: %f\n", mean(Kn_arr));

%% Get Kp and Kn from Groups 2 and 3  - Down/Up Trials

% Get subset of trials with left/right movement
moveDirNums = [trials.MoveDir];
md_ind = moveDirNums == 2;
trials_md = trials(md_ind);

% Get subset of trials with group number 0 or 1
groupNumbers = [trials_md.GroupNumber];
gn2 = groupNumbers == 2;
gn3 = groupNumbers == 3;
gn23 = gn2 | gn3;

trials_gn23 = trials_md(gn23);

% Calculate average Kp, Kn for trials
Kp_arr = [];
Kn_arr = [];
for i = 1:length(trials_gn23)
    [Kp, Kn] = GetKs(trials_gn23(i).Data.xdot, trials_gn23(i).Data.xdotdot);
    Kp_arr = [Kp_arr, Kp];
    Kn_arr = [Kn_arr, Kn];
end

fprintf("Groups 2 & 3\n");
fprintf("Kp: %f\n", mean(Kp_arr));
fprintf("Kn: %f\n", mean(Kn_arr));


%% Subject 1 Data Corrections

for i = 1:length(trials)
    if (trials(i).MoveDir == 2)         % Down - Up
        if (trials(i).TargetDirNum == 1)
            trials(i).TargetDirNum = 3;  % Left becomes down
        elseif (trials(i).TargetDirNum == 2)
            trials(i).TargetDirNum = 4; % Right becomes up
        end
    end
end


for i = 1:length(trials)
    if (trials(i).TargetDirNum == 3)
        trials(i).TargetDirText = 'Down';
    elseif (trials(i).TargetDirNum == 4)
        trials(i).TargetDirText = 'Up';
    end
end

%% EMG List 
emglist = [1,2,3,5,6,8];

for i = 1:length(emglist)
    emg = emglist(i);
    figure
    plot(trials(2).Data.Emg(emg,:))
end


