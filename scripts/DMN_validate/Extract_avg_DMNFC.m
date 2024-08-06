% extract DMN from baseline

% read in subj-session-dose correspondence
subSeshDose=readtable('~/subjSeshDoseCorresp.csv');
% add paths
ToolFolder='/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti/tool_folder';
addpath(genpath(ToolFolder));
% get subj list
subjPrefix=repmat('sub-MDMA0',17,1);
subjSuffix=["01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17"];
subjList=strcat(subjPrefix,subjSuffix')
% initialize cross-task FC
overallFC_left = [];
overallFC_right = [];
% for collecting each task
for task=["rs1" "rs2" "emotion" "gambling" "wm"]
	% set common fp
	commonFP=['/scratch/users/apines/data/mdma/'];
	% initialize output vectors: left and right each vertex
	outL=zeros(2562,17);
	outR=zeros(2562,17);
	% for each subj except 4 6 and 10
	for s=[1 2 3 5 7 8 9 11 12 13 14 15 16 17]
		disp(s)
	        % get session info
	        seshInfo=subSeshDose{s,2:5};
	        bvFP=[commonFP subjList(s) '/' seshInfo{1} '/' subjList(s) '_' seshInfo{1} '_' task '_DMNFC_map.mat'];
	        bvFP=strjoin(bvFP,'');
		% load in baseline csvs
		bv=load(bvFP).CorStruct;
		overallFC_left=[overallFC_left bv.DMNcor_L];
                overallFC_right=[overallFC_right bv.DMNcor_R];
	end
end
% combine across tasks
meanFC_left=mean(overallFC_left,2);
meanFC_right=mean(overallFC_right,2);
% set NaN's to 0 and print
meanFC_left(isnan(meanFC_left))=0;
meanFC_right(isnan(meanFC_right))=0;
Vis_FC(meanFC_left,meanFC_right,'~/MeanDMN_FC.png')
