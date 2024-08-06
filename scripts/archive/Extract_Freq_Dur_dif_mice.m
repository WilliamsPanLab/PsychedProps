function Extract_ang_dif_mice()
% extract network props from each session for this frequency band
% add paths
ToolFolder='/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti/tool_folder';
addpath(genpath(ToolFolder));
% set subj list, missing run 6 for 7507
subjList = {'m2000','m7507','m7520','m7522','m7589','m7594'};;
% initialize output file, 6 sessions, 1 for name
outDF=cell(length(subjList),7);
outDF_dur=cell(length(subjList),7);
% set common fp
commonFP=['/scratch/users/apines/data/mouse/'];
% for each mouse
for s=1:length(subjList);
	disp(s)
	for run=1:6
		% insert subject ID as first column
		outDF(s,1)=subjList(s);
		outDF_dur(s,1)=subjList(s);
	        ofFP_base='/scratch/users/apines/data/mouse/'
		ofFP=strjoin([ofFP_base subjList(s) '_' num2str(run) '_BU_counts.mat'],'');
		ofFP2=strjoin([ofFP_base subjList(s) '_' num2str(run) '_BU_meanDur.mat'],'');
		% load in baseline csvs
		if exist(ofFP,'file')
			file=load(ofFP);
	        	outDF(s,(run+1))={file.BUPcount};
			file2=load(ofFP2);
			outDF_dur(s,(run+1))={file2.meanBUPdur};
		end
	% end for each run
	end
% end for each mouse
end
% save out matrix
writecell(outDF,'/oak/stanford/groups/leanew1/users/apines/data/mice_AngFreq_Merged_mice_LSD.csv')
writecell(outDF_dur,'/oak/stanford/groups/leanew1/users/apines/data/mice_AngDur_Merged_mice_LSD.csv')
end
