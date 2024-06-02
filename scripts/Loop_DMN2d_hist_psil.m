% generate 2d dmn histograms (radial plots) for each subject session and task, to be merged later 
% list psil subjects
subjects={'PS03','PS16','PS18','PS19','PS21','PS24','PS93','PS96','PS98'}
subjects={'PS21','PS24','PS93','PS96','PS98'}
subjects={'PS96'}
% list all possible psil sessions (will just return a lot of "file is missing")
sessions = {'Baseline1','Baseline2','Baseline3','Baseline4','Baseline5','Baseline6','Baseline7','Baseline8','Baseline9',...
		'Drug1','Drug2'...
	       	'Between1','Between2','Between3','Between4','Between5',...
	       	'After1','After2','After3','After4','After5','After6','After7','After8'};
tasks = {'rs1', 'rs2','rs3','rs4','rs5','rs6'};
% for each subject
for i=1:length(subjects)
	subj = subjects{i}
	% for each session
	for j = 1:length(sessions)
		sesh = sessions{j};
		% set filepath that things would be saved to
		outFP=['/scratch/users/apines/data/psil/' subj '/' sesh];
		% and make an output directory that is needed
		system(['mkdir /oak/stanford/groups/leanew1/users/apines/data/p50/' subj]);
		% for each task
		for k = 1:length(tasks)
			task = tasks{k};
			% if file exists, then continue
			filename=[outFP '/' subj '_' sesh '_' task '_k1_Prop_Feats_gro.csv']
			if isfile(filename)
				% Interpolate Time series to in-between points and faces
				InterpolateTS_psil(subj,sesh,task)
				% Combine facewise time series with DMN angles into one .mat
				Combine_FacewiseTS_psil(subj,sesh,task)
				% make respective directory
				system(['mkdir /oak/stanford/groups/leanew1/users/apines/data/p50/' subj '/' sesh]);
				system(['mkdir /oak/stanford/groups/leanew1/users/apines/data/p50/' subj '/' sesh '/figs/']);
				% try using system to run python script to print out 2d histogram and radar visualization
				system(['python Viz_AngMag_psil.py ' subj ' ' sesh ' ' task]);
			else
				disp(['File is missing:' subj ' ' sesh ' ' task]);
			end
		end
	end
end
