% generate 2d dmn histograms (radial plots) for each subject session and task, to be merged later 
subjects = {'sub-MDMA001', 'sub-MDMA002', 'sub-MDMA003', 'sub-MDMA005', 'sub-MDMA007', 'sub-MDMA008', 'sub-MDMA009', 'sub-MDMA011', 'sub-MDMA012', 'sub-MDMA013', 'sub-MDMA014', 'sub-MDMA015', 'sub-MDMA016', 'sub-MDMA017'};
sessions = {'ses-00', 'ses-01', 'ses-02', 'ses-03'};
tasks = {'rs1', 'rs2', 'emotion', 'gambling', 'wm'};
% for each subject
for i=1:length(subjects)
	subj = subjects{i}
	% for each session
	for j = 1:length(sessions)
		sesh = sessions{j};
		% set filepath that things would be saved to
		outFP=['/scratch/users/apines/data/mdma/' subj '/' sesh];
		% for each task
		for k = 1:length(tasks)
			task = tasks{k};
			% if file exists, then continue
			filename=[outFP '/' subj '_' sesh '_' task '_k1_Prop_Feats_gro.csv']
			if isfile(filename)
				% Interpolate Time series to in-between points and faces
				InterpolateTS(subj,sesh,task)
				% Combine facewise time series with DMN angles into one .mat
				Combine_FacewiseTS(subj,sesh,task)
				% try using system to run python script to print out 2d histogram and radar visualization
				system(['python Viz_AngMag.py ' subj ' ' sesh ' ' task]);
			else
				disp(['File is missing:' subj ' ' sesh ' ' task]);
			end
		end
	end
end
