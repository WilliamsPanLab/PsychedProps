function Extract_ang_dif_mice_Diaz()
% extract network props from each session for this frequency band
% add paths
ToolFolder='/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti/tool_folder';
addpath(genpath(ToolFolder));
% set subj list
subjList={'m7531','m7576','m7520','m7522','m7590'};
% initialize output file, 6 sessions, 1 for name
outDF=cell(length(subjList),7);
% set common fp
commonFP=['/scratch/users/apines/data/mouse/'];
% for each mouse
for s=1:length(subjList);
	disp(s)
	for run=1:5
		% insert subject ID as first column
		outDF(s,1)=subjList(s);
	        ofFP='/scratch/users/apines/data/mouse/Diaz/'
		ofFP=strjoin([ofFP subjList(s) '_' num2str(run) '_Prop_Feats_Mag.csv'],'');
		% load in baseline csvs
		if exist(ofFP,'file')
			file=readmatrix(ofFP);
	        	outDF(s,(run+1))=num2cell(file(1,2));
		end
	% end for each run
	end
% end for each mouse
end
% save out matrix
writecell(outDF,'/oak/stanford/groups/leanew1/users/apines/data/mice_magsMerged_mice_Diaz.csv')
end
