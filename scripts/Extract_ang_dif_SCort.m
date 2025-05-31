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
	% initialize output table: difference by Subcort struct, for 80 and 120,
	outDF=zeros(17,20);
	% set common fp
	commonFP=['/scratch/users/apines/data/mdma/'];
	% for each subj except 4 6 10
	for s=[1 2 3 5 7 8 9 11 12 13 14 15 16 17]
		disp(s)
	        % get session info
	        seshInfo=subSeshDose{s,2:5};
	        bvFP_CL=[commonFP subjList(s) '/' seshInfo{1} '/' subjList(s) '_' seshInfo{1} '_' task '_CaudL_Prop_Feats_gro.csv'];
	        bvFP_CL=strjoin(bvFP_CL,'');
		bvFP_CR=[commonFP subjList(s) '/' seshInfo{1} '/' subjList(s) '_' seshInfo{1} '_' task '_CaudR_Prop_Feats_gro.csv'];
	        bvFP_CR=strjoin(bvFP_CR,'');
		bvFP_HL=[commonFP subjList(s) '/' seshInfo{1} '/' subjList(s) '_' seshInfo{1} '_' task '_HippoL_Prop_Feats_gro.csv'];
	        bvFP_HL=strjoin(bvFP_HL,'');
		bvFP_HR=[commonFP subjList(s) '/' seshInfo{1} '/' subjList(s) '_' seshInfo{1} '_' task '_HippoR_Prop_Feats_gro.csv'];
		bvFP_HR=strjoin(bvFP_HR,'');
		% and grab other sessions
		pFP_CL=[commonFP subjList(s) '/' seshInfo{2} '/' subjList(s) '_' seshInfo{2} '_' task '_CaudL_Prop_Feats_gro.csv'];
	        pFP_CL=strjoin(pFP_CL,'');
		pFP_CR=[commonFP subjList(s) '/' seshInfo{2} '/' subjList(s) '_' seshInfo{2} '_' task '_CaudR_Prop_Feats_gro.csv'];
		pFP_CR=strjoin(pFP_CR,'');
		pFP_HL=[commonFP subjList(s) '/' seshInfo{2} '/' subjList(s) '_' seshInfo{2} '_' task '_HippoL_Prop_Feats_gro.csv'];
	        pFP_HL=strjoin(pFP_HL,'');
		pFP_HR=[commonFP subjList(s) '/' seshInfo{2} '/' subjList(s) '_' seshInfo{2} '_' task '_HippoR_Prop_Feats_gro.csv'];
		pFP_HR=strjoin(pFP_HR,'');
		m1FP_CL=[commonFP subjList(s) '/' seshInfo{3} '/' subjList(s) '_' seshInfo{3} '_' task '_CaudL_Prop_Feats_gro.csv'];
	        m1FP_CL=strjoin(m1FP_CL,'');
		m1FP_CR=[commonFP subjList(s) '/' seshInfo{3} '/' subjList(s) '_' seshInfo{3} '_' task '_CaudR_Prop_Feats_gro.csv'];
		m1FP_CR=strjoin(m1FP_CR,'');
		m1FP_HL=[commonFP subjList(s) '/' seshInfo{3} '/' subjList(s) '_' seshInfo{3} '_' task '_HippoL_Prop_Feats_gro.csv'];
	        m1FP_HL=strjoin(m1FP_HL,'');
		m1FP_HR=[commonFP subjList(s) '/' seshInfo{3} '/' subjList(s) '_' seshInfo{3} '_' task '_HippoR_Prop_Feats_gro.csv'];
		m1FP_HR=strjoin(m1FP_HR,'');
		m2FP_CL=[commonFP subjList(s) '/' seshInfo{4} '/' subjList(s) '_' seshInfo{4} '_' task '_CaudL_Prop_Feats_gro.csv'];
	        m2FP_CL=strjoin(m2FP_CL,'');
		m2FP_CR=[commonFP subjList(s) '/' seshInfo{4} '/' subjList(s) '_' seshInfo{4} '_' task '_CaudR_Prop_Feats_gro.csv'];
		m2FP_CR=strjoin(m2FP_CR,'');
		m2FP_HL=[commonFP subjList(s) '/' seshInfo{4} '/' subjList(s) '_' seshInfo{4} '_' task '_HippoL_Prop_Feats_gro.csv'];
	        m2FP_HL=strjoin(m2FP_HL,'');
		m2FP_HR=[commonFP subjList(s) '/' seshInfo{4} '/' subjList(s) '_' seshInfo{4} '_' task '_HippoR_Prop_Feats_gro.csv'];
		m2FP_HR=strjoin(m2FP_HR,'');
		% load in baseline csvs
		if exist(bvFP_CL,'file')
		bvCL=readmatrix(bvFP_CL);
		bvCR=readmatrix(bvFP_CR);
		bvHL=readmatrix(bvFP_HL);
		bvHR=readmatrix(bvFP_HR);
	        outDF(s,1)=bvCL(1,2);
		outDF(s,2)=bvCR(1,2);
		outDF(s,3)=bvHL(1,2);
		outDF(s,4)=bvHR(1,2);
		% load in remaining frames
	        childfp=['/scratch/users/apines/data/mdma/' subjList(s) '/' seshInfo{1} ];
                bvCSIfp=[childfp '/' subjList(s) '_' seshInfo{1} '_task-' task '_ValidSegments_Trunc.txt'];
                bvCSIfp=strjoin(bvCSIfp,'');
                CSI = importdata(bvCSIfp);
                % trailing -1 is because the count column (,2) is inclusive of the start TR (,1)
                numTRsVS=sum(CSI(:,2));
                outDF(s,17)=numTRsVS;	
		end
		% load in placebo csvs
		if exist(pFP_CL,'file')
		pCL=readmatrix(pFP_CL);
		pCR=readmatrix(pFP_CR);
		pHL=readmatrix(pFP_HL);
		pHR=readmatrix(pFP_HR);
		% extract in subcort structs: placebo
		outDF(s,5)=pCL(1,2);
		outDF(s,6)=pCR(1,2);
		outDF(s,7)=pHL(1,2);
		outDF(s,8)=pHR(1,2);
		% load in remaining frames
		childfp=['/scratch/users/apines/data/mdma/' subjList(s) '/' seshInfo{2} ];
                pCSIfp=[childfp '/' subjList(s) '_' seshInfo{2} '_task-' task '_ValidSegments_Trunc.txt'];
                pCSIfp=strjoin(pCSIfp,'');
                CSI = importdata(pCSIfp);
                numTRsVS=sum(CSI(:,2));
                outDF(s,18)=numTRsVS;
		end
		% load mdma 1 (80 mg)
		if exist(m1FP_CL,'file')
		m1CL=readmatrix(m1FP_CL);
		m1CR=readmatrix(m1FP_CR);
		m1HL=readmatrix(m1FP_HL);
		m1HR=readmatrix(m1FP_HR);
		% extract in subcort structs: mdma 1
		outDF(s,9)=m1CL(1,2);
		outDF(s,10)=m1CR(1,2);
		outDF(s,11)=m1HL(1,2);
		outDF(s,12)=m1HR(1,2);
		% load in remaining frames
		childfp=['/scratch/users/apines/data/mdma/' subjList(s) '/' seshInfo{3} ];
                m1CSIfp=[childfp '/' subjList(s) '_' seshInfo{3} '_task-' task '_ValidSegments_Trunc.txt'];
                m1CSIfp=strjoin(m1CSIfp,'');
                CSI = importdata(m1CSIfp);
                numTRsVS=sum(CSI(:,2));
                outDF(s,19)=numTRsVS;
		end
		% load mdma 2 (120 mg)
		if exist(m2FP_CL,'file')
		m2CL=readmatrix(m2FP_CL);
		m2CR=readmatrix(m2FP_CR);
		m2HL=readmatrix(m2FP_HL);
		m2HR=readmatrix(m2FP_HR);
		% extract in subcort structs: mdma 2
		outDF(s,13)=m2CL(1,2);
		outDF(s,14)=m2CR(1,2);
		outDF(s,15)=m2HL(1,2);
		outDF(s,16)=m2HR(1,2);
		childfp=['/scratch/users/apines/data/mdma/' subjList(s) '/' seshInfo{4} ];
                m2CSIfp=[childfp '/' subjList(s) '_' seshInfo{4} '_task-' task '_ValidSegments_Trunc.txt'];
                m2CSIfp=strjoin(m2CSIfp,'');
                CSI = importdata(m2CSIfp);
                numTRsVS=sum(CSI(:,2));
                outDF(s,20)=numTRsVS;
		end
	end
	% save out matrix
	writematrix(outDF,strjoin(['/oak/stanford/groups/leanew1/users/apines/data/' task '_SCort_propsMerged.csv'],''))
end
