function Calc_AvgBup(subj)
% read in subj-session-dose correspondence
restoredefaultpath
subSeshDose=readtable('~/subjSeshDoseCorresp.csv');
% add paths
addpath(genpath('/oak/stanford/groups/leanew1/users/apines/libs/'))
addpath('/oak/stanford/groups/leanew1/users/apines/scripts/')
% Load in surface data
surfL = ['/oak/stanford/groups/leanew1/users/apines/surf/lh.sphere'];
surfR = ['/oak/stanford/groups/leanew1/users/apines/surf/rh.sphere'];
% surface topography
[vx_l, faces_l] = read_surf(surfL);
[vx_r, faces_r] = read_surf(surfR);
% +1 the faces: begins indexing at 0
faces_l = faces_l + 1;
faces_r = faces_r + 1;
% faces_L
F_L=faces_l;
% vertices V
V_L=vx_l;
% faces_R
F_R=faces_r;
% vertices V
V_R=vx_r;

% load in medial wall + SNR so we don't have to loop over every single vertex and then mask out later
% add TSNR mask, includes medial wall
mwAndTSNR_L='/oak/stanford/groups/leanew1/users/apines/fs4surf/lh.Mask_SNR.func.gii';
mwAndTSNR_R='/oak/stanford/groups/leanew1/users/apines/fs4surf/rh.Mask_SNR.func.gii';
mwAndTSNR_L=gifti(mwAndTSNR_L).cdata(:,1);
mwAndTSNR_R=gifti(mwAndTSNR_R).cdata(:,1);
mw_L=zeros(1,2562);
mw_L(mwAndTSNR_L==1)=1;
mw_R=zeros(1,2562);
mw_R(mwAndTSNR_R==1)=1;
% convert to faces
F_MW_L=sum(mw_L(faces_l),2)./3;
F_MW_R=sum(mw_R(faces_r),2)./3;
% convert "partial" medial wall to medial wall
F_MW_L=ceil(F_MW_L);
F_MW_R=ceil(F_MW_R);
% face mask indices
fmwIndVec_l=find(F_MW_L);
fmwIndVec_r=find(F_MW_R);
% load in DMN to make more thorough mask
networks=load(['/oak/stanford/groups/leanew1/users/apines/data/Atlas_Visualize/gro_Nets_fs4.mat']);
%% k = 1 to select DMN.
Dnet_LH=networks.nets.Lnets(:,1);
Dnet_RH=networks.nets.Rnets(:,1);
nets_LH=networks.nets.Lnets(:,1);
nets_RH=networks.nets.Rnets(:,1);
% create face-wise network mask
DMN_bool_L=sum(nets_LH(faces_l),2)./3;
DMN_bool_R=sum(nets_RH(faces_r),2)./3;
DMN_bool_L(DMN_bool_L>.3)=1;
DMN_bool_R(DMN_bool_R>.3)=1;
DMN_bool_L(DMN_bool_L<.3)=0;
DMN_bool_R(DMN_bool_R<.3)=0;
DMN_bool_L=logical(DMN_bool_L);
DMN_bool_R=logical(DMN_bool_R);
% combine with medial wall mask
MasterMask_L=DMN_bool_L;
MasterMask_R=DMN_bool_R;
MasterMask_L(fmwIndVec_l)=0;
MasterMask_R(fmwIndVec_r)=0;
% should be 1116 faces for left, 996 for right
mw_L=MasterMask_L;
mw_R=MasterMask_R;

% get subj list
subjPrefix=repmat('sub-MDMA0',17,1);
subjSuffix=["01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17"];
subjList=strcat(subjPrefix,subjSuffix')

% set common filepath
commonFP=['/scratch/users/apines/data/mdma/'];

% initialize cross-condition distance vectors
pl_80_L=zeros(sum(mw_L),1);
pl_80_R=zeros(sum(mw_R),1);
pl_120_L=zeros(sum(mw_L),1);
pl_120_R=zeros(sum(mw_R),1);
pl_drug_L=zeros(sum(mw_L),1);
pl_drug_R=zeros(sum(mw_R),1);
bv_80_L=zeros(sum(mw_L),1);
bv_80_R=zeros(sum(mw_R),1);
bv_120_L=zeros(sum(mw_L),1);
bv_120_R=zeros(sum(mw_R),1);
bv_drug_L=zeros(sum(mw_L),1);
bv_drug_R=zeros(sum(mw_R),1);
bv_plac_L=zeros(sum(mw_L),1);
bv_plac_R=zeros(sum(mw_R),1);
m80_120_L=zeros(sum(mw_L),1);
m80_120_R=zeros(sum(mw_R),1);

% to inherit same code as group-level scripting
s=find(strcmp(subj,subjList));
% now set subj to acutal subject name for filepaths
% should be redundant with input
subj=subjList(s);
% get session-condition correspondence
seshInfo=subSeshDose{s,2:5};
% set conditions
conditions={'bv','pl','m1','m2'};
% initialize opfl structure to hold each condition/hemisphere
OpFl_rs1 = struct();
OpFl_rs2 = struct();

% for rs1
task='rs1'
% now loop over each condition to load in each and concatenate resting-state angular time series
for c=1:4
	condition=conditions{c};
	% LOAD IN ANGULAR TIME SERIES instead (and facewise)
	fpl=[commonFP '/' subj '/' seshInfo{c} '/' subj '_' seshInfo{c} '_' task '_k1_Prop_TS_dmn_L.csv'];
	fpr=[commonFP '/' subj '/' seshInfo{c} '/' subj '_' seshInfo{c} '_' task '_k1_Prop_TS_dmn_R.csv'];
	% need to add # or surviving TRs check to EVERY load
        survivingTrsFP=strjoin([commonFP '/' subj '/' seshInfo{c} '/' subj '_' seshInfo{c} '_task-' task '_ValidSegments_Trunc.txt'],'');
        survivingTrs=load(survivingTrsFP);
        % first condition needed because some desolate scans dont even have a second column for # of valid TRs, because no sequences of TRs meeting motion criteria exists
	if size(survivingTrs, 2) > 1 && sum(survivingTrs(:,2))>250;
		OpFl_rs1.(condition).L=dlmread(strjoin(fpl,''));
		OpFl_rs1.(condition).R=dlmread(strjoin(fpr,''));
	else
		OpFl_rs1.(condition).L=0;
		OpFl_rs1.(condition).R=0;
	end
end

% for rs2
task='rs2'
% now loop over each condition to load in each and concatenate resting-state angular time series
for c=1:4
        condition=conditions{c};
	fpl=[commonFP '/' subj '/' seshInfo{c} '/' subj '_' seshInfo{c} '_' task '_k1_Prop_TS_dmn_L.csv'];
	fpr=[commonFP '/' subj '/' seshInfo{c} '/' subj '_' seshInfo{c} '_' task '_k1_Prop_TS_dmn_R.csv'];
        % need to add # or surviving TRs check to EVERY load
        survivingTrsFP=strjoin([commonFP '/' subj '/' seshInfo{c} '/' subj '_' seshInfo{c} '_task-' task '_ValidSegments_Trunc.txt'],'');
        survivingTrs=load(survivingTrsFP);
        % first condition needed because some desolate scans dont even have a second column for # of valid TRs, because no sequences of TRs meeting motion criteria exists
        if size(survivingTrs, 2) > 1 && sum(survivingTrs(:,2))>250;
                OpFl_rs2.(condition).L=dlmread(strjoin(fpl,''));
                OpFl_rs2.(condition).R=dlmread(strjoin(fpr,''));
        else    
                OpFl_rs2.(condition).L=0;
                OpFl_rs2.(condition).R=0;
        end     
end 

% now concatenate the scans into a master OpFl struct
OpFl=struct();
for c=1:4
	% for this condition
	condition=conditions{c};	
	% if rs1 and rs2 are ~=0, concatenate
	if OpFl_rs1.(condition).L(1) ~=0 && OpFl_rs2.(condition).L(1) ~=0
		OpFl.(condition).L=cat(2, OpFl_rs1.(condition).L, OpFl_rs2.(condition).L);
	% if rs2 is 0, but not rs1, use rs1
	elseif OpFl_rs1.(condition).L(1) ~=0 && OpFl_rs2.(condition).L(1) ==0
		OpFl.(condition).L=OpFl_rs1.(condition).L;
	% if rs1 is 0 but not rs2, use rs2
	elseif OpFl_rs1.(condition).L(1)==0 && OpFl_rs2.(condition).L(1) ~=0
		OpFl.(condition).L=OpFl_rs2.(condition).L;
	end
	% right hemisphere
	% if rs1 and rs2 are ~=0, concatenate
        if OpFl_rs1.(condition).R(1) ~=0 && OpFl_rs2.(condition).R(1) ~=0
                OpFl.(condition).R=cat(2, OpFl_rs1.(condition).R, OpFl_rs2.(condition).R);
        % if rs2 is 0, but not rs1, use rs1
        elseif OpFl_rs1.(condition).R(1) ~=0 && OpFl_rs2.(condition).R(1) ==0
                OpFl.(condition).R=OpFl_rs1.(condition).R;
        % if rs1 is 0 but not rs2, use rs2
        elseif OpFl_rs1.(condition).R(1) ==0 && OpFl_rs2.(condition).R(1) ~=0
		OpFl.(condition).R=OpFl_rs2.(condition).R;
	end
            
end

% saveout filepath
outFP=['/scratch/users/apines/data/mdma/' subj];

% for each face, get average BUP angle instead of time series of them
if OpFl.bv.L(1)>0
	% BV
	bv_Angles_L=OpFl.bv.L;
	bv_Angles_R=OpFl.bv.R;
	% convert to % BUP
	BV_L=(sum(bv_Angles_L < 90,2) / size(bv_Angles_L,2) * 100);
	BV_R=(sum(bv_Angles_R < 90,2) / size(bv_Angles_R,2) * 100);
	save(strjoin([outFP '/AvgBup_BV_L.mat'],""),'BV_L');
	save(strjoin([outFP '/AvgBup_BV_R.mat'],""),'BV_R');
end
% placebo
if OpFl.pl.L(1)>0
	% PL
	pl_Angles_L=OpFl.pl.L;
	pl_Angles_R=OpFl.pl.R;
	PL_L=(sum(pl_Angles_L < 90,2) / size(pl_Angles_L,2) * 100);
	PL_R=(sum(pl_Angles_R < 90,2) / size(pl_Angles_R,2) * 100);
	save(strjoin([outFP '/AvgBup_PL_L.mat'],""),'PL_L');
	save(strjoin([outFP '/AvgBup_PL_R.mat'],""),'PL_R');
end
% m80
if OpFl.m1.L(1)>0
	% m1
	m1_Angles_L=OpFl.m1.L;
	m1_Angles_R=OpFl.m1.R;
	M1_L=(sum(m1_Angles_L < 90,2) / size(m1_Angles_L,2) * 100);
	M1_R=(sum(m1_Angles_R < 90,2) / size(m1_Angles_R,2) * 100);
	% save out average BUP to scratch
        save(strjoin([outFP '/AvgBup_M1_L.mat'],""),'M1_L');
        save(strjoin([outFP '/AvgBup_M1_R.mat'],""),'M1_R');
end
% m120
if OpFl.m2.L(1)>0
	% m2
	m2_Angles_L=OpFl.m2.L;
	m2_Angles_R=OpFl.m2.R;
	M2_L=(sum(m2_Angles_L < 90,2) / size(m2_Angles_L,2) * 100);
	M2_R=(sum(m2_Angles_R < 90,2) / size(m2_Angles_R,2) * 100);
	% save out average BUP to scratch
	save(strjoin([outFP '/AvgBup_M2_L.mat'],""),'M2_L');
	save(strjoin([outFP '/AvgBup_M2_R.mat'],""),'M2_R');
end

