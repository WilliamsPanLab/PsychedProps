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

% label filepath common for each subj
commonFP=['/scratch/users/apines/data/mdma/'];

% get subj list
subjPrefix=repmat('sub-MDMA0',17,1);
subjSuffix=["01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17"];
subjList=strcat(subjPrefix,subjSuffix')

% for each subj except 4
for s=[1 2 3 5 6 7 8 9 10 11 12 13 14 15 16 17]
	% get session info
	seshInfo=subSeshDose{s,2:5};
	bvFP=[commonFP subjList(s) '/' seshInfo{1} '/' subjList(s) '_' seshInfo{1} '_PGGDist_rs.mat'];
	bvFP=strjoin(bvFP,'')
	pFP=[commonFP subjList(s) '/' seshInfo{2} '/' subjList(s) '_' seshInfo{2} '_PGGDist_rs.mat'];
	pFP=strjoin(pFP,'')
	m1FP=[commonFP subjList(s) '/' seshInfo{3} '/' subjList(s) '_' seshInfo{3} '_PGGDist_rs.mat'];
	m1FP=strjoin(m1FP,'')
	m2FP=[commonFP subjList(s) '/' seshInfo{4} '/' subjList(s) '_' seshInfo{4} '_PGGDist_rs.mat'];
	m2FP=strjoin(m2FP,'')
	% load in baseline session
	bv=load(bvFP);
	% load in placebo session
	p1=load(pFP);
	% mdma sessions
	m1=load(m1FP);
	m2=load(m2FP);
	% concat sober into one ts
	sobLeft=vertcat(bv.AngDist.gLeft,p1.AngDist.gLeft);
	sobRight=vertcat(bv.AngDist.gRight,p1.AngDist.gRight);
	% concat mdma into one ts
	mLeft=vertcat(m1.AngDist.gLeft,m2.AngDist.gLeft);
	mRight=vertcat(m1.AngDist.gRight,m2.AngDist.gRight);
	% init tvec and pvec
	tvec_L=zeros(length(g_noMW_combined_L),1);
	pvec_L=zeros(length(g_noMW_combined_L),1);
	tvec_R=zeros(length(g_noMW_combined_R),1);
	pvec_R=zeros(length(g_noMW_combined_R),1);
	counter=0;
	% for each left face
	for F=g_noMW_combined_L
		counter=counter+1;
		% t-test
		[h,p,ci,stats] = ttest(sobLeft(:,F),mLeft(:,F));
		tvec_L(counter)=stats.tstat;
		pvec_L(counter)=p;
	end
	counter=0;
	% for each right face
	for F=g_noMW_combined_R
		counter=counter+1;
		% t-test
                [h,p,ci,stats] = ttest(sobRight(:,F),mRight(:,F));
                tvec_R(counter)=stats.tstat;
                pvec_R(counter)=p;
	end
	% combine ps
	ps=[pvec_L' pvec_R'];
	% combine ts
	ts=[tvec_L' tvec_R'];
	% mafdr ps
	fdred=mafdr(ps);
	% use it to thresh T's
	ts(fdred>0.05)=0;
	ts_L=ts(1:length(g_noMW_combined_L));
	ts_R=ts((length(g_noMW_combined_L)+1):length(ts));
	% png out filename
	pngFN=['~/' subjList(s) '_ts.png'];
	% visfacevec of threshed T's
	Vis_Facevec(ts_L,ts_R,strjoin(pngFN,''))
end

