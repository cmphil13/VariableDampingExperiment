function files = GetAllFilesInDir(base_dir)
%GETALLFILESINDIR returns an array of all files in the directory inputted
%as well as all files in all subdirectories

% Get all directories and files in base_dir
dirinfo = dir(base_dir);
% Get rid of '.' and '..'
for i = length(dirinfo):-1:1
    if (strcmp(dirinfo(i).name,'.') || strcmp(dirinfo(i).name, '..'))
        dirinfo(i) = [];
    end
end

% Get files in directory
fileInds = ([dirinfo.isdir] == 0);
files = dirinfo(fileInds);

% Recursively search subdirectories
dirInds = [dirinfo.isdir] == 1;
nDirs = sum(dirInds);

if (nDirs ~= 0)
    subdirs = dirinfo(dirInds);
    for i = 1:nDirs
        subdir = fullfile(subdirs(i).folder, subdirs(i).name);
        files = [files; GetAllFilesInDir(subdir)];
    end
end

end

