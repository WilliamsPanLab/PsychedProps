% read in subj-session-dose correspondence
subSeshDose=readtable('~/subjSeshDoseCorresp.csv');

% load in fixed mw mask
addpath(genpath('/oak/stanford/groups/leanew1/users/apines/libs'))

%%% Load in surface data
SubjectsFolder = '/oak/stanford/groups/leanew1/users/apines/fs5surf';
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
mw_L=zeros(1,10242);
mw_L(mwIndVec_l)=1;
mw_R=zeros(1,10242);
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
g_noMW_combined_L=setdiff([1:20480],fmwIndVec_l);
g_noMW_combined_R=setdiff([1:20480],fmwIndVec_r);
noMW_L=setdiff(1:10242,mwIndVec_l);
noMW_R=setdiff(1:10242,mwIndVec_r);

% label filepath common for each subj
commonFP=['/scratch/users/apines/data/mdma/'];

% get subj list
subjPrefix=repmat('sub-MDMA0',17,1);
subjSuffix=["01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17"];
subjList=strcat(subjPrefix,subjSuffix')

% make a mean T's vector for analyses after this script
meanTs=zeros(2,17);
% for each subj except 4
for s=[1 2 3 5 6 7 8 9 10 11 12 13 14 15 16 17]
	% get session info
	seshInfo=subSeshDose{s,2:5};
	bvFP=[commonFP subjList(s) '/' seshInfo{1} '/OpFl_timeseries_L_fs5.mat'];
	bvFP=strjoin(bvFP,'')
	bvFPr=[commonFP subjList(s) '/' seshInfo{1} '/OpFl_timeseries_R_fs5.mat'];
        bvFPr=strjoin(bvFPr,'')
	pFP=[commonFP subjList(s) '/' seshInfo{2} '/OpFl_timeseries_L_fs5.mat'];
	pFP=strjoin(pFP,'')
	pFPr=[commonFP subjList(s) '/' seshInfo{2} '/OpFl_timeseries_R_fs5.mat'];
        pFPr=strjoin(pFPr,'')
	m1FP=[commonFP subjList(s) '/' seshInfo{3} '/OpFl_timeseries_L_fs5.mat'];
	m1FP=strjoin(m1FP,'')
	m1FPr=[commonFP subjList(s) '/' seshInfo{3} '/OpFl_timeseries_R_fs5.mat'];
        m1FPr=strjoin(m1FPr,'')
	m2FP=[commonFP subjList(s) '/' seshInfo{4} '/OpFl_timeseries_L_fs5.mat'];
	m2FP=strjoin(m2FP,'')
	m2FPr=[commonFP subjList(s) '/' seshInfo{4} '/OpFl_timeseries_R_fs5.mat'];
        m2FPr=strjoin(m2FPr,'')
	% load in baseline session
	bv=load(bvFP);
	% pull out time series of interest
	bv=squeeze(bv.LvertTS_num(:,3,:));
	bvr=load(bvFPr);
	bvr=squeeze(bvr.RvertTS_num(:,3,:));
	% load in placebo session
	p1=load(pFP);
	p1=squeeze(p1.LvertTS_num(:,3,:));
	p1r=load(pFPr);
	p1r=squeeze(p1r.RvertTS_num(:,3,:));
	% mdma sessions
	m1=load(m1FP);
	m1=squeeze(m1.LvertTS_num(:,3,:));
	m1r=load(m1FPr);
	m1r=squeeze(m1r.RvertTS_num(:,3,:));
	m2=load(m2FP);
	m2=squeeze(m2.LvertTS_num(:,3,:));
	m2r=load(m2FPr);
	m2r=squeeze(m2r.RvertTS_num(:,3,:));
	% concat sober into one ts
	sobLeft=horzcat(bv,p1);
	sobRight=horzcat(bvr,p1r);
	%sobLeft=p1;
	%sobRight=p1r;
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
		[h,p,ci,stats] = ttest(sobLeft(F,:),mLeft(F,:));
		tvec_L(counter)=stats.tstat;
		pvec_L(counter)=p;
	end
	counter=0;
	% for each right face
	for F=noMW_R
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
	% mafdr ps
	fdred=mafdr(ps);
	% print pre-thresh average of ts, positive reflects further distance from PG in unintox.
	meanTs(1,s)=s
	meanTs(2,s)=mean(ts(~isnan(ts)))
	% use it to thresh T's
	ts(fdred>0.05)=0;
	ts_L=ts(1:length(noMW_L));
	ts_R=ts((length(noMW_L)+1):length(ts));
	% png out filename
	pngFN=['~/' subjList(s) '_ts.png'];
	% visfacevec of threshed T's
	Vis_VertvecFs5(ts_L,ts_R,strjoin(pngFN,''))
end
% saveout mean ts
csvwrite('~/fs5_ts.csv',meanTs)
