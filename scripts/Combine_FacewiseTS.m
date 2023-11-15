function Combine_FacewiseTS(subj,sesh,task)
% for table reading
restoredefaultpath
% load in interpolated facewise time series, combine with angular distance time series
childfp=['/scratch/users/apines/data/mdma/' subj '/' sesh ];
ofpl=[childfp '/' subj '_' sesh '_task-' task '_p2mm_masked_interp_L_faces_DMN.csv'];
ofpr=[childfp '/' subj '_' sesh '_task-' task '_p2mm_masked_interp_R_faces_DMN.csv'];
% load in interpTS
iTS_l=readtable(ofpl);
iTS_r=readtable(ofpr);
% read in angular distances
aTS_l=readtable([childfp '/' subj '_' sesh '_' task '_Prop_TS_dmn_L.csv']);
aTS_r=readtable([childfp '/' subj '_' sesh '_' task '_Prop_TS_dmn_R.csv']);
% for freesurfer functions
addpath(genpath('/oak/stanford/groups/leanew1/users/apines/libs/'))
% get dimensions to combine on
spatialDim=size(iTS_l);
spatialDim=spatialDim(1);
temporalDim=size(iTS_l);
temporalDim=temporalDim(2);
% combine
FullMatrix_L=zeros(spatialDim,temporalDim,2);
FullMatrix_L(:,:,1)=table2array(iTS_l);
FullMatrix_L(:,:,2)=table2array(aTS_l);
% right hemi
spatialDim=size(iTS_r);
spatialDim=spatialDim(1);
temporalDim=size(iTS_r);
temporalDim=temporalDim(2);
FullMatrix_R=zeros(spatialDim,temporalDim,2);
FullMatrix_R(:,:,1)=table2array(iTS_r);
FullMatrix_R(:,:,2)=table2array(aTS_r);
% Load in surface data
%addpath(genpath('/oak/stanford/groups/leanew1/users/apines/libs'))
%SubjectsFolder = '/oak/stanford/groups/leanew1/users/apines/surf';
%surfL = [SubjectsFolder '/lh.sphere'];
%surfR = [SubjectsFolder '/rh.sphere'];
% surface topography
%[vx_l, faces_l] = read_surf(surfL);
%[vx_r, faces_r] = read_surf(surfR);
% +1 the faces: begins indexing at 0
%faces_l = faces_l + 1;
%faces_r = faces_r + 1;
% faces_L
%F_L=faces_l;
% vertices V
%V_L=vx_l;
% faces_R
%F_R=faces_r;
% vertices V
%V_R=vx_r;
% use native freesurfer command for mw mask indices
%surfML = [SubjectsFolder '/lh.Medial_wall.label'];
%mwIndVec_l = read_medial_wall_label(surfML);
%surfMR = [SubjectsFolder '/rh.Medial_wall.label'];
%mwIndVec_r = read_medial_wall_label(surfMR);
% make binary "is medial wall" vector for vertices
%mw_L=zeros(1,2562);
%mw_L(mwIndVec_l)=1;
%mw_R=zeros(1,2562);
%mw_R(mwIndVec_r)=1;
% convert to faces
%F_MW_L=sum(mw_L(faces_l),2)./3;
%F_MW_R=sum(mw_R(faces_r),2)./3;
% convert "partial" medial wall to medial wall
%F_MW_L=ceil(F_MW_L);
%F_MW_R=ceil(F_MW_R);
% face mask indices
%fmwIndVec_l=find(F_MW_L);
%fmwIndVec_r=find(F_MW_R);
%g_noMW_combined_L=setdiff([1:5120],fmwIndVec_l);
%g_noMW_combined_R=setdiff([1:5120],fmwIndVec_r);
% load in InclLeft InclRight
%Incl=load('/oak/stanford/groups/leanew1/users/apines/fs5surf/medial_wall_nullGrad_vectors.mat');
%InclLeft=Incl.InclLeft;
%InclRight=Incl.InclRight;
% load in DMN
%networks=load(['/oak/stanford/groups/leanew1/users/apines/data/RobustInitialization/group_Nets_fs4.mat']);
%nets_LH=networks.nets.Lnets;
%nets_RH=networks.nets.Rnets;
%DMN_L=nets_LH(:,2);
%DMN_R=nets_RH(:,2);
% convert to faces
%DMN_L_f=sum(DMN_L(faces_l),2)./3;
%DMN_R_f=sum(DMN_R(faces_r),2)./3;
% mask DMN
%DMN_L_f=DMN_L_f(g_noMW_combined_L);
%DMN_R_f=DMN_R_f(g_noMW_combined_R);
%DMN_L_f=DMN_L_f(InclLeft);
%DMN_R_f=DMN_R_f(InclRight);
% now use masked DMN to make a DMN mask with equivalent faces to time series
%DMN_L_f_mask=DMN_L_f>.3;
%DMN_R_f_mask=DMN_R_f>.3;
% now maks the time series
%FullMatrix_L=FullMatrix_L(DMN_L_f_mask,:,:);
%FullMatrix_R=FullMatrix_R(DMN_R_f_mask,:,:);
% saveout
ofl=[childfp '/' subj '_' sesh '_task-' task '_p2mm_masked_Bold_and_Angles_L.mat'];
ofr=[childfp '/' subj '_' sesh '_task-' task '_p2mm_masked_Bold_and_Angles_R.mat'];
save(ofl, 'FullMatrix_L');
save(ofr, 'FullMatrix_R');
