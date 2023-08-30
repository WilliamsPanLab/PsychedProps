function InterpolateTS(subj,sesh)
% this function is to interpolate the functional time series to within-segment and between-TR timepoints. Should be v stable in band-passed fMR signal
addpath(genpath('/oak/stanford/groups/leanew1/users/apines/libs/'))
% load in time series (downsampled)
childfp=['/scratch/users/apines/data/mdma/' subj '/' sesh ];
% load in data
fpL=[childfp '/' subj '_' sesh '_task-rs_p2mm_masked_L.mgh'];
fpR=[childfp '/' subj '_' sesh '_task-rs_p2mm_masked_R.mgh'];
dataL=MRIread(fpL);
dataR=MRIread(fpR);
% squeeze to get rid of extra dimensions
TRs_l_g=squeeze(dataL.vol);
TRs_r_g=squeeze(dataR.vol);
% get sizes
sizefMR=size(TRs_l_g);
lengthTS=sizefMR(2);
% load in continuous segment information
CSIfp=[childfp '/' subj '_' sesh '_task-rs_ValidSegments_Trunc.txt'];
CSI = importdata(CSIfp);
% assure that TR count is the same between time series and valid segments txt
SegNum=size(CSI);
SegNum=SegNum(1);
% trailing -1 is because the count column (,2) is inclusive of the start TR (,1)
numTRsVS=CSI(SegNum,1)+CSI(SegNum,2)-1;
if numTRsVS ~= lengthTS
        disp('TRs from Valid Segments txt and mgh file do not match. Fix it.')
        return
end
% store original timepoints as sequence
originalTimepoints = 1:lengthTS; 
%%% derive segment-ending TRs from CSI to factor them out
% get the number that each segment starts with
SEtrs=CSI(:,1);
% take out first starting point (i.e., first TR) from this exclusion list. Only seeking to omit interpolation between book-ends of segments
SEtrs=SEtrs(2:end);
% and finally, -1 to get 1 prior to starting point
SEtrs=SEtrs-1;
% now omit these TRs from originalTimepoints in terms of what TRs we want to interpolate between
% AKA remove segment-ending TRs from this sequence to not interpolate between-segment TRs
timepointsToInterpolate = setdiff(originalTimepoints, SEtrs);
% Calculate halfway timepoints
halfwayTimepoints = timepointsToInterpolate(1:end-1) + 0.5;
% initialize interpolated data
InterpData_L=zeros(sizefMR(1),length(halfwayTimepoints));
InterpData_R=zeros(sizefMR(1),length(halfwayTimepoints));
% Interpolate data at halfway timepoints
for v=1:10242
	interpolatedData_L(v,:) = interp1(originalTimepoints, TRs_l_g(v,:), halfwayTimepoints);
	interpolatedData_R(v,:) = interp1(originalTimepoints, TRs_r_g(v,:), halfwayTimepoints);
end
% save it out
ofpl=[childfp '/' subj '_' sesh '_task-rs_p2mm_masked_interp_L.mgh'];
ofpr=[childfp '/' subj '_' sesh '_task-rs_p2mm_masked_interp_R.mgh'];
% saveout
dataL.vol=interpolatedData_L;
dataR.vol=interpolatedData_R;
MRIwrite(dataL,ofpl)
MRIwrite(dataR,ofpr)
% now convert vertices to faces: will have to save out as .csv outside of .mgh framework
SubjectsFolder = '/oak/stanford/groups/leanew1/users/apines/fs5surf';
% for surface data
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
F_R=faces_r;
% vertices V
V_L=vx_l;
% vertices V
V_R=vx_r;
% get incenters of triangles
TR_L = TriRep(F_L,V_L);
P_L = TR_L.incenters;
TR_R = TriRep(F_R,V_R);
P_R = TR_R.incenters;
% initialize output face matrix
tsLength=size(interpolatedData_L);
tsLength=tsLength(2);
faceOutL=zeros(20480,tsLength);
faceOutR=zeros(20480,tsLength);
% convert to faces
for t=1:tsLength
	interpolatedData_L_frame=interpolatedData_L(:,t);
	faceOutL(:,t)=sum(interpolatedData_L_frame(faces_l),2)./3;
	% right hemisphere
	interpolatedData_R_frame=interpolatedData_R(:,t);
	faceOutR(:,t)=sum(interpolatedData_R_frame(faces_r),2)./3;
end
% medial wall mask
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% use native freesurfer command for mw mask indices
surfML = '/oak/stanford/groups/leanew1/users/apines/fs5surf/lh.Medial_wall.label';
mwIndVec_l = read_medial_wall_label(surfML);
surfMR = '/oak/stanford/groups/leanew1/users/apines/fs5surf/rh.Medial_wall.label';
mwIndVec_r = read_medial_wall_label(surfMR);
% make binary "is medial wall" vector for vertices
mw_L=zeros(1,10242);
mw_L(mwIndVec_l)=1;
mw_R=zeros(1,10242);
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
g_noMW_combined_L=setdiff([1:20480],fmwIndVec_l);
g_noMW_combined_R=setdiff([1:20480],fmwIndVec_r);
%%%%%% mask it
faceOutL=faceOutL(g_noMW_combined_L,:);
faceOutR=faceOutR(g_noMW_combined_R,:);
%%%%% equivalent masking to extract_relativeangle.m
networks=load(['/oak/stanford/groups/leanew1/users/apines/data/RobustInitialization/group_Nets_fs5.mat']);
nets_LH=networks.nets.Lnets;
nets_RH=networks.nets.Rnets;
% network of interest
n_LH=nets_LH(:,2);
n_RH=nets_RH(:,2);
% calculate network gradients on sphere
ng_L = grad(F_L, V_L, n_LH);
ng_R = grad(F_R, V_R, n_RH);
% use medial wall mask as common starting point (from which to mask both opfl vecs and net grads further)
ng_L=ng_L(g_noMW_combined_L,:);
ng_R=ng_R(g_noMW_combined_R,:);
sumLeft=sum(ng_L,2);
sumRight=sum(ng_R,2);
% finds 0s in left and right network gradients
emptyLeft=find(~sumLeft);
emptyRight=find(~sumRight);
InclLeft=find(sumLeft);
InclRight=find(sumRight);
% last mask with inclleft and incl right on interp timeseries
faceOutL=faceOutL(InclLeft,:);
faceOutR=faceOutR(InclRight,:);
% saveout
ofpl=[childfp '/' subj '_' sesh '_task-rs_p2mm_masked_interp_L_faces.csv'];
ofpr=[childfp '/' subj '_' sesh '_task-rs_p2mm_masked_interp_R_faces.csv'];
% saveout
csvwrite(ofpl,faceOutL)
csvwrite(ofpr,faceOutR)
