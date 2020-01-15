function res = LoadEmgDataFile(emgFile)
    if (isfile(emgFile))
        [filepath, name, ext] = fileparts(emgFile);
        if (strcmp(ext, '.txt'))
            try
                res = load(emgFile);
            catch
                file = emgFile;
                % delete the last line of the text file
                [path, file1, ext] = fileparts(file);
                file1 = strcat(file1, '_2', ext);
                file1 = fullfile(path, file1);

                % copy emg file to EmgData_2.txt and delete original
                % EmgData.txt
                copyfile(file, file1);
                delete(file);

                % Get line count of file
                fid1 = fopen(file1);

                tline = fgetl(fid1);
                lineCount = 1;
                while ischar(tline)
                    tline = fgetl(fid1);
                    lineCount = lineCount + 1;
                end

                fclose(fid1);

                % copy all but last line to new EmgData.txt
                fid = fopen(file, 'w');
                fid1 = fopen(file1);
                for i = 1:lineCount - 2
                    tline = fgetl(fid1);
                    fprintf(fid, '%s\n', tline);
                end

                fclose(fid);
                fclose(fid1);
                
                % load fixed file
                res = load(file);
            end
                
        elseif (stcmp(ext, '.h5'))
            res = h5read(emgFile, '/EMG');
        end
    else
        res = [];
    end
end