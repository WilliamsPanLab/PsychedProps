% extract network props from each session

% add paths
ToolFolder='/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti/tool_folder';
addpath(genpath(ToolFolder));
% read in subj-session-dose correspondence
subSeshDose=readtable('~/subjSeshDoseCorresp.csv');
% get subj list
subjPrefix=repmat('sub-MDMA0',17,1);
subjSuffix=["01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17"];
subjList=strcat(subjPrefix,subjSuffix')
% for collecting each task
for task=["rs1" "rs2" "emotion" "gambling" "wm"]
	% initialize output table: difference in alff by each network, for 80 and 120, and in group and individualized atlases
	outDF=zeros(17,8);
	% set common fp
	commonFP=['/scratch/users/apines/data/mdma/'];
	% for each subj except 4 6 10
	for s=[1 2 3 5 7 8 9 11 12 13 14 15 16 17]
		disp(s)
	        % get session info NOTE K4 HAS ALL 4 VALUES, K1 has INTENDED VALUE AS WELL
	        seshInfo=subSeshDose{s,2:5};
	        bvFP=[commonFP subjList(s) '/' seshInfo{1} '/' subjList(s) '_' seshInfo{1} '_' task '_k4_Prop_Feats_gro.csv'];
	        bvFP=strjoin(bvFP,'');
		% and grab reverse phase encoding direction nifti
		pFP=[commonFP subjList(s) '/' seshInfo{2} '/' subjList(s) '_' seshInfo{2} '_' task '_k4_Prop_Feats_gro.csv'];
	        pFP=strjoin(pFP,'');
	        m1FP=[commonFP subjList(s) '/' seshInfo{3} '/' subjList(s) '_' seshInfo{3} '_' task '_k4_Prop_Feats_gro.csv'];
	        m1FP=strjoin(m1FP,'');
	        m2FP=[commonFP subjList(s) '/' seshInfo{4} '/' subjList(s) '_' seshInfo{4} '_' task '_k4_Prop_Feats_gro.csv'];
	        m2FP=strjoin(m2FP,'');
		% load in baseline csvs
		if exist(bvFP,'file')
		% check using props, but load in Dseg
		bvFP=[commonFP subjList(s) '/' seshInfo{1} '/' subjList(s) '_' seshInfo{1} '_' task '_DMNMag_TD.csv'];
		bvFP=strjoin(bvFP,'');
		bv=readmatrix(bvFP);
	        outDF(s,1)=bv(1,2);
		% load in remaining frames
	        childfp=['/scratch/users/apines/data/mdma/' subjList(s) '/' seshInfo{1} ];
                bvCSIfp=[childfp '/' subjList(s) '_' seshInfo{1} '_task-' task '_ValidSegments_Trunc.txt'];
                bvCSIfp=strjoin(bvCSIfp,'');
                CSI = importdata(bvCSIfp);
                % trailing -1 is because the count column (,2) is inclusive of the start TR (,1)
                numTRsVS=sum(CSI(:,2));
                outDF(s,5)=numTRsVS;	
		end
		% load in placebo csvs
		if exist(pFP,'file')
		% check using props, but load in Dseg
                pFP=[commonFP subjList(s) '/' seshInfo{2} '/' subjList(s) '_' seshInfo{2} '_' task '_DMNMag_TD.csv'];
                pFP=strjoin(pFP,'');
		p=readmatrix(pFP);
		% extract in group consensus atlas: placebo
		outDF(s,2)=p(1,2);
		childfp=['/scratch/users/apines/data/mdma/' subjList(s) '/' seshInfo{2} ];
                pCSIfp=[childfp '/' subjList(s) '_' seshInfo{2} '_task-' task '_ValidSegments_Trunc.txt'];
                pCSIfp=strjoin(pCSIfp,'');
                CSI = importdata(pCSIfp);
                numTRsVS=sum(CSI(:,2));
                outDF(s,6)=numTRsVS;
		end
		% load mdma 1 (80 mg)
		if exist(m1FP,'file')
		% check using props, but load in Dseg
                m1FP=[commonFP subjList(s) '/' seshInfo{3} '/' subjList(s) '_' seshInfo{3} '_' task '_DMNMag_TD.csv'];
                m1FP=strjoin(m1FP,'');
		m1=readmatrix(m1FP);
		% extact in group consensus atlas: mdma 1
		outDF(s,3)=m1(1,2);
		childfp=['/scratch/users/apines/data/mdma/' subjList(s) '/' seshInfo{3} ];
                m1CSIfp=[childfp '/' subjList(s) '_' seshInfo{3} '_task-' task '_ValidSegments_Trunc.txt'];
                m1CSIfp=strjoin(m1CSIfp,'');
                CSI = importdata(m1CSIfp);
                numTRsVS=sum(CSI(:,2));
                outDF(s,7)=numTRsVS;
		end
		% load mdma 2 (120 mg)
		if exist(m2FP,'file')
		% check using props, but load in Dseg
		m2FP=[commonFP subjList(s) '/' seshInfo{4} '/' subjList(s) '_' seshInfo{4} '_' task '_DMNMag_TD.csv'];
                m2FP=strjoin(m2FP,'');
		m2=readmatrix(m2FP);
		% extract in group consensus atlas: mdma 2
		outDF(s,4)=m2(1,2);
		childfp=['/scratch/users/apines/data/mdma/' subjList(s) '/' seshInfo{4} ];
                m2CSIfp=[childfp '/' subjList(s) '_' seshInfo{4} '_task-' task '_ValidSegments_Trunc.txt'];
                m2CSIfp=strjoin(m2CSIfp,'');
                CSI = importdata(m2CSIfp);
                numTRsVS=sum(CSI(:,2));
                outDF(s,8)=numTRsVS;
		end
	end
	% save out matrix
	writematrix(outDF,strjoin(['/oak/stanford/groups/leanew1/users/apines/data/' task '_DMNMag_TD_Merged.csv'],''))
end
