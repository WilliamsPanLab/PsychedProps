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
	% concat sober into one ts
	sobLeft=horzcat(bv,p1);
	sobRight=horzcat(bvr,p1r);
	% concat mdma into one ts
	mLeft=horzcat(m1,m2);
	mRight=horzcat(m1r,m2r);

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%%%%%%%% LOAD IN THAL FC DATA



	% no need for T's! we'll get two big vectors of delta FC and delta props

	% init tvec and pvec
	tvec_L=zeros(length(g_noMW_combined_L),1);
	pvec_L=zeros(length(g_noMW_combined_L),1);
	tvec_R=zeros(length(g_noMW_combined_R),1);
	pvec_R=zeros(length(g_noMW_combined_R),1);
	counter=0;
	
	%% note: will need to add DMN angles together for each face touching a vertex and to aggregate a vertex-wide measurement

	% for each left vert
	for F=1:length(g_noMW_combined_L)
		counter=counter+1;
		% t-test
		[h,p,ci,stats] = ttest(sobLeft(F,:),mLeft(F,:));
		tvec_L(counter)=stats.tstat;
		pvec_L(counter)=p;
	end
	counter=0;
	% for each right face
	for F=1:length(g_noMW_combined_R)
		counter=counter+1;
		% t-test
                [h,p,ci,stats] = ttest(sobRight(F,:),mRight(F,:));
                tvec_R(counter)=stats.tstat;
                pvec_R(counter)=p;
	end
	% combine ps
	ps=[pvec_L' pvec_R'];
	% combine ts
	ts=[tvec_L' tvec_R'];
	%%%%%%% mask out dmn
	% load in DMN
	networks=load(['/oak/stanford/groups/leanew1/users/apines/data/RobustInitialization/group_Nets_fs4.mat']);
	nets_LH=networks.nets.Lnets(:,2);
	nets_RH=networks.nets.Rnets(:,2);
	% convert to faces
	F_MW_L=sum(nets_LH(faces_l),2)./3;
	F_MW_R=sum(nets_RH(faces_r),2)./3;

	%%% should probably lessen this threshold to allow thalamic FC maps to fit in better

	% boolean out areas where .3 or less DMN
	F_MW_L(F_MW_L<0.3)=0;
	F_MW_L=ceil(F_MW_L);
	F_MW_R(F_MW_R<0.3)=0;
	F_MW_R=ceil(F_MW_R);
	% convert to within medial wall mask
	DMN_L=logical(F_MW_L(g_noMW_combined_L));
	DMN_R=logical(F_MW_R(g_noMW_combined_R));
	% mask ts and ps
	tvec_L_sub=tvec_L(DMN_L);
	tvec_R_sub=tvec_R(DMN_R);
	pvec_L_sub=pvec_L(DMN_L);
	pvec_R_sub=pvec_R(DMN_R);
	% combine ps
	ps=[pvec_L_sub' pvec_R_sub'];
	% combine ts
	ts=[tvec_L_sub' tvec_R_sub'];
	% mafdr ps
	fdred=mafdr(ps);
	% print pre-thresh average of ts, positive reflects further distance from DMNG in unintox.
	mean(ts)
	% use it to thresh T's
	ts(fdred>0.05)=0;
	% data needs to be mw-mask-length for each
	tvec_fin_l=zeros(1,length(g_noMW_combined_L));
	tvec_fin_r=zeros(1,length(g_noMW_combined_R));
	tvec_fin_l(DMN_L)=ts(1:sum(DMN_L));
	tvec_fin_r(DMN_R)=ts((sum(DMN_L)+1):(sum(DMN_L)+sum(DMN_R)));
	%%%%
	% png out filename
	pngFN=['~/' subjList(s) '_ts.png'];
	% visfacevec of threshed T's
	Vis_FaceVec(tvec_fin_l,tvec_fin_r,strjoin(pngFN,''))
end

