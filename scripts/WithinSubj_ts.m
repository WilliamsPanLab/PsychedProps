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
noMW_L=setdiff(1:2562,mwIndVec_l);
noMW_R=setdiff(1:2562,mwIndVec_r);

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
	bvFP=[commonFP subjList(s) '/' seshInfo{1} '/OpFl_timeseries_L.mat'];
	bvFP=strjoin(bvFP,'')
	bvFPr=[commonFP subjList(s) '/' seshInfo{1} '/OpFl_timeseries_R.mat'];
        bvFPr=strjoin(bvFPr,'')
	pFP=[commonFP subjList(s) '/' seshInfo{2} '/OpFl_timeseries_L.mat'];
	pFP=strjoin(pFP,'')
	pFPr=[commonFP subjList(s) '/' seshInfo{2} '/OpFl_timeseries_R.mat'];
        pFPr=strjoin(pFPr,'')
	m1FP=[commonFP subjList(s) '/' seshInfo{3} '/OpFl_timeseries_L.mat'];
	m1FP=strjoin(m1FP,'')
	m1FPr=[commonFP subjList(s) '/' seshInfo{3} '/OpFl_timeseries_R.mat'];
        m1FPr=strjoin(m1FPr,'')
	m2FP=[commonFP subjList(s) '/' seshInfo{4} '/OpFl_timeseries_L.mat'];
	m2FP=strjoin(m2FP,'')
	m2FPr=[commonFP subjList(s) '/' seshInfo{4} '/OpFl_timeseries_R.mat'];
        m2FPr=strjoin(m2FPr,'')
	% load in baseline session
	bv=load(bvFP);
	% pull out time series of interest
	bv=squeeze(bv.LvertTS(:,2,:));
	bvr=load(bvFPr);
	bvr=squeeze(bvr.RvertTS(:,2,:));
	% load in placebo session
	p1=load(pFP);
	p1=squeeze(p1.LvertTS(:,2,:));
	p1r=load(pFPr);
	p1r=squeeze(p1r.RvertTS(:,2,:));
	% mdma sessions
	m1=load(m1FP);
	m1=squeeze(m1.LvertTS(:,2,:));
	m1r=load(m1FPr);
	m1r=squeeze(m1r.RvertTS(:,2,:));
	m2=load(m2FP);
	m2=squeeze(m2.LvertTS(:,2,:));
	m2r=load(m2FPr);
	m2r=squeeze(m2r.RvertTS(:,2,:));
	% concat sober into one ts
	sobLeft=horzcat(bv,p1);
	sobRight=horzcat(bvr,p1r);
	% concat mdma into one ts
	mLeft=horzcat(m1,m2);
	mRight=horzcat(m1r,m2r);
	% init tvec and pvec
	tvec_L=zeros(length(noMW_L),1);
	pvec_L=zeros(length(noMW_L),1);
	tvec_R=zeros(length(noMW_R),1);
	pvec_R=zeros(length(noMW_R),1);
	counter=0;
	% for each left vert
	for F=noMW_L
		counter=counter+1;
		% t-test
		[h,p,ci,stats] = ttest(cell2mat(sobLeft(F,:)),cell2mat(mLeft(F,:)));
		tvec_L(counter)=stats.tstat;
		pvec_L(counter)=p;
	end
	counter=0;
	% for each right face
	for F=noMW_R
		counter=counter+1;
		% t-test
                [h,p,ci,stats] = ttest(cell2mat(sobRight(F,:)),cell2mat(mRight(F,:)));
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
	ts_L=ts(1:length(noMW_L));
	ts_R=ts((length(noMW_L)+1):length(ts));
	% png out filename
	pngFN=['~/' subjList(s) '_ts.png'];
	% visfacevec of threshed T's
	Vis_Vertvec(ts_L,ts_R,strjoin(pngFN,''))
end

