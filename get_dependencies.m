% 1. Get dependencies
deps = matlab.codetools.requiredFilesAndProducts('BATCH_HCP.m');

% 2. Define destination folder (change to your preferred path)
destFolder = fullfile(pwd, 'BATCH_HCP_Script_Dependencies');
if ~exist(destFolder, 'dir')
    mkdir(destFolder);
end

% 3. Copy all resolvable files

    filesList = cellstr(deps); % Ensure consistent indexing
    copiedCount = 0;
    
    for i = 1:length(filesList)
        srcFile = filesList{i};
        
        % Skip if file doesn't actually exist on disk
        if ~exist(srcFile, 'file')
            fprintf('⚠️  Skipped (not found): %s\n', srcFile);
            continue;
        end
        
        % Extract just the filename
        [~, name, ext] = fileparts(srcFile);
        destFile = fullfile(destFolder, [name ext]);
        
        % Handle naming collisions to avoid silent overwrites
        if exist(destFile, 'file')
            counter = 1;
            while exist(fullfile(destFolder, sprintf('%s_%d%s', name, counter, ext)), 'file')
                counter = counter + 1;
            end
            destFile = fullfile(destFolder, sprintf('%s_%d%s', name, counter, ext));
        end
        
        % Copy file
        copyfile(srcFile, destFile);
        copiedCount = copiedCount + 1;
    end
    
    fprintf('✅ Successfully copied %d files to: %s\n', copiedCount, destFolder);
