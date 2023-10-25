subjects = {'sub-MDMA001', 'sub-MDMA002', 'sub-MDMA003', 'sub-MDMA005', 'sub-MDMA007', 'sub-MDMA008', 'sub-MDMA009', 'sub-MDMA011', 'sub-MDMA012', 'sub-MDMA013', 'sub-MDMA014', 'sub-MDMA015', 'sub-MDMA016', 'sub-MDMA017'};
sessions = {'ses-00', 'ses-01', 'ses-02', 'ses-03'};

% Loop through subjects
for subjectCell = subjects
    % Loop through sessions
    for sessionCell = sessions
        % Extract the char from the cell
        subject = subjectCell{1};
        session = sessionCell{1}; 
	% Call the function with Extract_RelativeAngles(subject, session, 'rs1')
	% if output file doesn't exist
	if exist(['/scratch/users/apines/gp/PropFeats/' subject '_' session '_rs1_faceMatrix.mat'],'file')==0;
	Extract_RelativeAngles(subject, session, 'rs1');
	else
	end
	if exist(['/scratch/users/apines/gp/PropFeats/' subject '_' session '_rs2_faceMatrix.mat'],'file')==0;
	Extract_RelativeAngles(subject, session, 'rs2');
	else
	end
	if exist(['/scratch/users/apines/gp/PropFeats/' subject '_' session '_gambling_faceMatrix.mat'],'file')==0;
	Extract_RelativeAngles(subject, session, 'gambling');
	else
	end
	if exist(['/scratch/users/apines/gp/PropFeats/' subject '_' session '_emotion_faceMatrix.mat'],'file')==0;
	Extract_RelativeAngles(subject, session, 'emotion');
	else
	end
	if exist(['/scratch/users/apines/gp/PropFeats/' subject '_' session '_wm_faceMatrix.mat'],'file')==0;
	Extract_RelativeAngles(subject, session, 'wm');
	else
	end
    end
end
