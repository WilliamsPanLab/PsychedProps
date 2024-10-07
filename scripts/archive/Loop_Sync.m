% generate 2d dmn histograms (radial plots) for each subject session and task, to be merged later 
subjects = {'sub-MDMA001', 'sub-MDMA002', 'sub-MDMA003', 'sub-MDMA005', 'sub-MDMA007', 'sub-MDMA008', 'sub-MDMA009', 'sub-MDMA011', 'sub-MDMA012', 'sub-MDMA013', 'sub-MDMA014', 'sub-MDMA015', 'sub-MDMA016', 'sub-MDMA017'};
% done with first 7, 8 needs xcp re-run because scratch deleted it (temporary line of code if that isn't intuitive)
%subjects = {'sub-MDMA009', 'sub-MDMA011', 'sub-MDMA012', 'sub-MDMA013', 'sub-MDMA014', 'sub-MDMA015', 'sub-MDMA016', 'sub-MDMA017'};
% running 8 sep. because it just finished re-running (scratch mechanism deleted it)
%subjects={'sub-MDMA008'};
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
			% if file exists (optical flow ran), then continue
			filename=[outFP '/' subj '_' sesh '_' task '_k1_Prop_Feats_gro.csv']
			if isfile(filename)
				% convert optical flow output to vertices
				OpFl_toVerts(subj,sesh,task)
				% get intraparcel synchrony (left)
				IntraParcelSynchrony_ot_left(subj,sesh,task)
			end
		end
	end
end
