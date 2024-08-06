% extract DMN from baseline

% add paths
ToolFolder='/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti/tool_folder';
addpath(genpath(ToolFolder));
% get subj list
subjPrefix=repmat('PS',11,1);
subjSuffix=["02","03","16","18","19","21","24","93","96","98","99"];
subjList=strcat(subjPrefix,subjSuffix')
% initialize cross-task FC
overallFC_left = [];
overallFC_right = [];
% for collecting each task
for task=["rs1" "rs2" "rs3" "rs4" "rs5" "rs6"]
	% set common fp
	commonFP=['/scratch/users/apines/data/psil/'];
	% initialize output vectors: left and right each vertex
	outL=zeros(2562,11,8);
	outR=zeros(2562,11,8);
	% for each subj except 2
        for s=[2 3 4 5 6 7 8 9 10 11];	
		% for up to 8 iterations
	        for i=1:8
			disp(s)
	        	% pull in baseline
			bvFP=[commonFP subjList(s) '/Baseline' num2str(i) '/' subjList(s) '_Baseline' num2str(i) '_' task '_DMNFC_map.mat'];
	        	bvFP=strjoin(bvFP,'');
			if exist(bvFP,'file')
				% load in baseline fcs
				bv=load(bvFP).CorStruct;
				overallFC_left=[overallFC_left bv.DMNcor_L];
				overallFC_right=[overallFC_right bv.DMNcor_R];
			end
		end
	end
end
% combine across tasks
meanFC_left=mean(overallFC_left,2);
meanFC_right=mean(overallFC_right,2);
% set NaN's to 0 and print
meanFC_left(isnan(meanFC_left))=0;
meanFC_right(isnan(meanFC_right))=0;
Vis_FC(meanFC_left,meanFC_right,'~/MeanDMN_FC_psil.png')
