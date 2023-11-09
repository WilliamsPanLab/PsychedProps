restoredefaultpath
% read in subj-session-dose correspondence
subSeshDose=readtable('~/subjSeshDoseCorresp.csv');

% load in fixed mw mask
addpath(genpath('/oak/stanford/groups/leanew1/users/apines/libs'))

%%% Load in surface data
SubjectsFolder = '/oak/stanford/groups/leanew1/users/apines/surf';
surfL = [SubjectsFolder '/lh.sphere'];
surfR = [SubjectsFolder '/rh.sphere'];
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
% use native freesurfer command for mw mask indices
surfML = [SubjectsFolder '/lh.Medial_wall.label'];
mwIndVec_l = read_medial_wall_label(surfML);
surfMR = [SubjectsFolder '/rh.Medial_wall.label'];
mwIndVec_r = read_medial_wall_label(surfMR);
% make binary "is medial wall" vector for vertices
mw_L=zeros(1,2562);
mw_L(mwIndVec_l)=1;
mw_R=zeros(1,2562);
mw_R(mwIndVec_r)=1;
% convert to faces
F_MW_L=sum(mw_L(faces_l),2)./3;
F_MW_R=sum(mw_R(faces_r),2)./3;
% convert "partial" medial wall to medial wall
F_MW_L=ceil(F_MW_L);
F_MW_R=ceil(F_MW_R);
% face mask indices
fmwIndVec_l=find(F_MW_L);
fmwIndVec_r=find(F_MW_R);
% make medial wall vector
g_noMW_combined_L=setdiff([1:5120],fmwIndVec_l);
g_noMW_combined_R=setdiff([1:5120],fmwIndVec_r);
noMW_L=setdiff(1:2562,mwIndVec_l);
noMW_R=setdiff(1:2562,mwIndVec_r);

% label filepath common for each subj
commonFP=['/oak/stanford/groups/leanew1/users/apines/data/gp/PropFeats/'];
FCfp=['/scratch/groups/leanew1/xcpd_outP50_36p_bp/xcp_d/'];

% get subj list
subjPrefix=repmat('sub-MDMA0',17,1);
subjSuffix=["01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17"];
subjList=strcat(subjPrefix,subjSuffix')

% make a mean T's vector for analyses after this script
% for each subj except 4
for s=[1 2 3 5 6 7 8 9 10 11 12 13 14 15 16 17]
	% get session info
	seshInfo=subSeshDose{s,2:5};
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%%%%% LOAD IN PROPAGATIONS DATA
	bvFP=[commonFP subjList(s) '/' subjList(s) '_' seshInfo{1} '_Prop_TS_dmn_L.csv'];
	bvFP=strjoin(bvFP,'')
	bvFPr=[commonFP subjList(s) '/' subjList(s) '_' seshInfo{1} '_Prop_TS_dmn_R.csv'];
        bvFPr=strjoin(bvFPr,'')
	pFP=[commonFP subjList(s) '/' subjList(s) '_' seshInfo{2} '_Prop_TS_dmn_L.csv'];
	pFP=strjoin(pFP,'')
	pFPr=[commonFP subjList(s) '/' subjList(s) '_' seshInfo{2} '_Prop_TS_dmn_R.csv'];
        pFPr=strjoin(pFPr,'')
	m1FP=[commonFP subjList(s) '/' subjList(s) '_' seshInfo{3} '_Prop_TS_dmn_L.csv'];
	m1FP=strjoin(m1FP,'')
	m1FPr=[commonFP subjList(s) '/' subjList(s) '_' seshInfo{3} '_Prop_TS_dmn_R.csv'];
        m1FPr=strjoin(m1FPr,'')
	m2FP=[commonFP subjList(s) '/' subjList(s) '_' seshInfo{4} '_Prop_TS_dmn_L.csv'];
	m2FP=strjoin(m2FP,'')
	m2FPr=[commonFP subjList(s) '/' subjList(s) '_' seshInfo{4} '_Prop_TS_dmn_R.csv'];
        m2FPr=strjoin(m2FPr,'')
	% load in baseline session
	bv=load(bvFP);
	bvr=load(bvFPr);
	% load in placebo session
	p1=load(pFP);
	p1r=load(pFPr);
	% mdma sessions
	m1=load(m1FP);
	m1r=load(m1FPr);
	m2=load(m2FP);
	m2r=load(m2FPr);
	% placebo into sober ts
	sobLeft=p1;
	sobRight=p1r;
	% concat mdma into one ts
	mLeft=horzcat(m1,m2);
	mRight=horzcat(m1r,m2r);
	% convert face-wise values to vertices by summing distributions at the faces each vertex touches per vertex
	% get length of time series
	tl=size(bv,2);
	tl_m=size(mLeft,2);
	% initialize out t-stats
	outt_L=zeros(1,2562);
	outt_R=zeros(1,2562);
	% for each vertex
	for v=1:2562
		% get faces that touch vertex
		faces_l=find(any(F_L==v,2));
		faces_r=find(any(F_R==v,2));
		% reconstruct time series with medial wall areas included
		reconTS_l=zeros(5120,tl);
		reconTS_r=zeros(5120,tl);
		reconTS_l(g_noMW_combined_L,:)=sobLeft;
		reconTS_r(g_noMW_combined_R,:)=sobRight;
		% get values at those faces
		sob_vals_l=reconTS_l(faces_l,:);
		sob_vals_r=reconTS_r(faces_r,:);
		% load in drug values
		reconTS_l=zeros(5120,tl_m);
		reconTS_r=zeros(5120,tl_m);
		reconTS_l(g_noMW_combined_L,:)=mLeft;
		reconTS_r(g_noMW_combined_R,:)=mRight;
		% get values at those faces
		m_vals_l=reconTS_l(faces_l,:);
		m_vals_r=reconTS_r(faces_r,:);
		% retrieve t-statistic for each vertex
		[h,p,ci,stats] = ttest2(sob_vals_l(:),m_vals_l(:));
		outt_L(v)=stats.tstat;
		[h,p,ci,stats] = ttest2(sob_vals_r(:),m_vals_r(:));
		outt_R(v)=stats.tstat;
		% remember to re-mask out medial wall!!
	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%%%%%%%% LOAD IN ALFF DATA
	bvFP=[FCfp subjList(s) '/' seshInfo{1} '/func/' subjList(s) '_' seshInfo{1} '_alff1_L_3k.func.gii'];
	bvFP=char(strjoin(bvFP,''))
	bvFPr=[FCfp subjList(s) '/' seshInfo{1} '/func/' subjList(s) '_' seshInfo{1} '_alff1_R_3k.func.gii'];
	bvFPr=char(strjoin(bvFPr,''))
	pFP=[FCfp subjList(s) '/' seshInfo{2} '/func/' subjList(s) '_' seshInfo{2} '_alff1_L_3k.func.gii'];
	pFP=char(strjoin(pFP,''))
	pFPr=[FCfp subjList(s) '/' seshInfo{2} '/func/' subjList(s) '_' seshInfo{2} '_alff1_R_3k.func.gii'];
	pFPr=char(strjoin(pFPr,''))
	m1FP=[FCfp subjList(s) '/' seshInfo{3} '/func/' subjList(s) '_' seshInfo{3} '_alff1_L_3k.func.gii'];
	m1FP=char(strjoin(m1FP,''))
	m1FPr=[FCfp subjList(s) '/' seshInfo{3} '/func/' subjList(s) '_' seshInfo{3} '_alff1_R_3k.func.gii'];
	m1FPr=char(strjoin(m1FPr,''))
	m2FP=[FCfp subjList(s) '/' seshInfo{4} '/func/' subjList(s) '_' seshInfo{4} '_alff1_L_3k.func.gii'];
	m2FP=char(strjoin(m2FP,''))
	m2FPr=[FCfp subjList(s) '/' seshInfo{4} '/func/' subjList(s) '_' seshInfo{4} '_alff1_R_3k.func.gii'];
	m2FPr=char(strjoin(m2FPr,''))
	bvFP2=[FCfp subjList(s) '/' seshInfo{1} '/func/' subjList(s) '_' seshInfo{1} '_alff2_L_3k.func.gii'];
        bvFP2=char(strjoin(bvFP,''))
        bvFPr2=[FCfp subjList(s) '/' seshInfo{1} '/func/' subjList(s) '_' seshInfo{1} '_alff2_R_3k.func.gii'];
        bvFPr2=char(strjoin(bvFPr,''))
        pFP2=[FCfp subjList(s) '/' seshInfo{2} '/func/' subjList(s) '_' seshInfo{2} '_alff2_L_3k.func.gii'];
        pFP2=char(strjoin(pFP,''))
        pFPr2=[FCfp subjList(s) '/' seshInfo{2} '/func/' subjList(s) '_' seshInfo{2} '_alff2_R_3k.func.gii'];
        pFPr2=char(strjoin(pFPr,''))
        m1FP2=[FCfp subjList(s) '/' seshInfo{3} '/func/' subjList(s) '_' seshInfo{3} '_alff2_L_3k.func.gii'];
        m1FP2=char(strjoin(m1FP,''))
        m1FPr2=[FCfp subjList(s) '/' seshInfo{3} '/func/' subjList(s) '_' seshInfo{3} '_alff2_R_3k.func.gii'];
        m1FPr2=char(strjoin(m1FPr,''))
        m2FP2=[FCfp subjList(s) '/' seshInfo{4} '/func/' subjList(s) '_' seshInfo{4} '_alff2_L_3k.func.gii'];
        m2FP2=char(strjoin(m2FP,''))
        m2FPr2=[FCfp subjList(s) '/' seshInfo{4} '/func/' subjList(s) '_' seshInfo{4} '_alff2_R_3k.func.gii'];
        m2FPr2=char(strjoin(m2FPr,''))
	% load in baseline session
	bvalf=gifti(bvFP).cdata;
	bvalf_r=gifti(bvFPr).cdata;
	bvalf2=gifti(bvFP2).cdata;
	bvalf2_r=gifti(bvFPr2).cdata;
	% load in placebo session
	p1alf=gifti(pFP).cdata;
	p1alf_r=gifti(pFPr).cdata;
	p1alf2=gifti(pFP2).cdata;
	p1alf2_r=gifti(pFPr2).cdata;
	% mdma sessions
	m1alf=gifti(m1FP).cdata;
	m1alf_r=gifti(m1FPr).cdata;
	m1alf2=gifti(m1FP2).cdata;
	m1alf2_r=gifti(m1FPr2).cdata;
	m2alf=gifti(m2FP).cdata;
	m2alf_r=gifti(m2FPr).cdata;
	m2alf2=gifti(m2FP2).cdata;
	m2alf2_r=gifti(m2FPr2).cdata;
	% average alff 1 and 2 together
	bvAlf=(bvalf+bvalf2)/2;
	bvAlf_r=(bvalf_r+bvalf2_r)/2;
	p1Alf=(p1alf+p1alf2)/2;
	p1Alf_r=(p1alf_r+p1alf2_r)/2;
	m1Alf=(m1alf+m1alf2)/2;
	m1Alf_r=(m1alf_r+m1alf2_r)/2;
	m2Alf=(m2alf+m2alf2)/2;
	m2Alf_r=(m2alf_r+m2alf2_r)/2;
	% placebo into as sober
	sobLeft_alff=p1Alf;
	sobRight_alff=p1Alf_r;
	% average mdma sessions
	mLeft_alff=(m1Alf+m2Alf)/2;
	mRight_alff=(m1Alf_r+m2Alf_r)/2;
	
	% difference FC vector
	diffAlff_L=sobLeft_alff-mLeft_alff;
	diffAlff_R=sobRight_alff-mRight_alff;

	%%%%%%% mask out dmn
	% load in DMN
	networks=load(['/oak/stanford/groups/leanew1/users/apines/data/RobustInitialization/group_Nets_fs4.mat']);
	nets_LH=networks.nets.Lnets(:,2);
	nets_RH=networks.nets.Rnets(:,2);
	% mask vertices
	nets_LH(nets_LH<0.1)=0;
	nets_RH(nets_RH<0.1)=0;
	% convert to within DMN mask
	ts_L=outt_L(logical(nets_LH));
	FCdelta_L=diffFC_L(logical(nets_LH))';
	ts_R=outt_R(logical(nets_RH));
	FCdelta_R=diffFC_R(logical(nets_RH))';
	% get nan indices to mask out
	nanBool=isnan(ts_L);
	ts_L=ts_L(~nanBool);
	FCdelta_L=FCdelta_L(~nanBool);
	nanBool=isnan(ts_R);
	ts_R=ts_R(~nanBool);
	FCdelta_R=FCdelta_R(~nanBool);
	% CREATE SCATTERPLOTS LEFT AND RIGHT	
	figure
	scatter(ts_L,FCdelta_L,'MarkerEdgeAlpha', 0.5)
	% png out filename
	pngFN=['~/' subjList(s) '_L_delta_fMR.png'];
	print(strjoin(pngFN,''),'-dpng')
	figure
        scatter(ts_R,FCdelta_R,'MarkerEdgeAlpha', 0.5)
        % png out filename
        pngFN=['~/' subjList(s) '_R_delta_fMR.png'];
        print(strjoin(pngFN,''),'-dpng')
	corr(ts_R',FCdelta_R')
end

