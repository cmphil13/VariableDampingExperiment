function res = LoadKukaDataFile(file)
    % Load data
    data = load(file);
    dataStruct = struct;
    
    % Split into functional parts
    dataStruct.JointAngle = data(:,1:7);
    dataStruct.Force = data(:,8:9);
    dataStruct.Pose = data(:,10:15);
    dataStruct.TargetDirection = data(:,16);
    dataStruct.Damping = data(:,17);
    dataStruct.xdot = data(:,18);
    dataStruct.xdotdot = data(:,19);
    if (size(data,2) > 19)
        dataStruct.MotorTorque = data(:,20:26);
    end
    
    res = dataStruct;
end