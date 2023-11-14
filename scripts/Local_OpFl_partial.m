% Define subjects and sessions
subjects = {'sub-MDMA015','sub-MDMA016', 'sub-MDMA017'};
sessions = {'ses-00', 'ses-01', 'ses-02', 'ses-03'};

% Loop over subjects
for subjIdx = 1:length(subjects)
    subj = subjects{subjIdx};
    
    % Loop over sessions
    for seshIdx = 1:length(sessions)
        sesh = sessions{seshIdx};
        
        fprintf('Processing subject %s, session %s\n', subj, sesh);

	% input dir
	inputDir = ['/scratch/users/apines/data/mdma/' subj '/' sesh];
        % output dir
	outputDir = ['/oak/stanford/groups/leanew1/users/apines/OpFlAngDs/mdma/' subj];
        mkdir(outputDir);

        % Tasks to process
        tasks = {'rs1', 'rs2', 'emotion', 'gambling', 'wm'};

        % Loop over tasks
        for taskIdx = 1:length(tasks)
            task = tasks{taskIdx};
            
            % File path to check
            filePath = sprintf('%s/%s_%s_%s_OpFl.mat', inputDir, subj, sesh, task)

            % Check if the file exists
            if ~exist(filePath, 'file')
                fprintf('File %s does not exist. Skipping...\n', filePath);
            else
                % Call your MATLAB function if the file doesn't exist
                Extract_RelativeAngles(subj, sesh, task);
            end
        end

        % Optionally, add a pause between function calls to avoid overloading the system
        pause(1);
    end
end

disp('OpFl complete');

