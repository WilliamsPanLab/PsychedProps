% Define subjects and sessions
subjects = {'sub-MDMA016', 'sub-MDMA017'};
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

        Extract_AmygFC(subj,sesh)

        % Optionally, add a pause between function calls to avoid overloading the system
        pause(1);
    end
end

disp('OpFl complete');

