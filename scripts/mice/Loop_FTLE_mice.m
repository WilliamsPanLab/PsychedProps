% generate 2d dmn histograms (radial plots) for each subject session and task, to be merged later 
%subjects = {'sub-MDMA001', 'sub-MDMA002', 'sub-MDMA003', 'sub-MDMA005', 'sub-MDMA007', 'sub-MDMA008', 'sub-MDMA009', 'sub-MDMA011', 'sub-MDMA012', 'sub-MDMA013', 'sub-MDMA014', 'sub-MDMA015', 'sub-MDMA016', 'sub-MDMA017'};
% done with first 7, 8 needs xcp re-run because scratch deleted it (temporary line of code if that isn't intuitive)
%subjects = {'sub-MDMA009', 'sub-MDMA011', 'sub-MDMA012', 'sub-MDMA013', 'sub-MDMA014', 'sub-MDMA015', 'sub-MDMA016', 'sub-MDMA017'};
% running 8 sep. because it just finished re-running (scratch mechanism deleted it)
%subjects={'m2000','m7507','m7520','m7589','m7594'};
subjects={'m7520'};
sessions = {1,2,3,4,5,6};
% for each subject
for i=1:length(subjects)
	subj = subjects{i}
	% for each session
	for j = 1:length(sessions)
		sesh = sessions{j};
		% set filepath that things would be saved to
		outFP=['/scratch/users/apines/data/mouse'];
		% if file exists, then continue
		filename=[outFP '/' subj '_' num2str(sesh) '_Prop_Feats_gro.csv']
		if isfile(filename)
			run_FlowPortrait_mice(subj,sesh)
		else
			disp(['File is missing:' subj ' ' num2str(sesh)]);
		end
	end
end
