function us = OpFl_mdma(subj,sesh,task)
% add paths
addpath(genpath('/oak/stanford/groups/leanew1/users/apines/libs/'))
%%%%%%%%%%%%%%%%%%%% Set parameters
N = 10; % Vector spherical harmonics degree
h = 1; % finite-difference height vector of triangles
alpha = 1; % scales Tikhonov regularization, > alpha = > spatiotemporal smoothness
s = 1; % R(u), regularizing functional, scales Tikhonov regularization more rapidly via penalization of first spatial and first temporal derivative
%%%%%%%%%%%%%%%%%%%%

childfp=['/scratch/users/apines/data/mdma/' subj '/' sesh ];

% load in data
fpL=[childfp '/' subj '_' sesh '_L_' task '_p2mm_masked_L.mgh'];
fpR=[childfp '/' subj '_' sesh '_R_' task '_p2mm_masked_R.mgh'];

dataL=MRIread(fpL).vol;
dataR=MRIread(fpR).vol;
% squeeze to get rid of extra dimensions
TRs_l_g=squeeze(dataL);
TRs_r_g=squeeze(dataR);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% load in continuous segment indices
CSIfp=[childfp '/' subj '_' sesh '_task-' task '_ValidSegments_Trunc.txt'];
CSI = importdata(CSIfp);
% assure that TR count is the same between time series and valid segments txt
SegNum=size(CSI);
SegNum=SegNum(1);
% number of trs in fmri ts
mr_ts_trs=size(TRs_l_g);
mr_ts_trs=mr_ts_trs(2);
% trailing -1 is because the count column (,2) is inclusive of the start TR (,1)
numTRsVS=CSI(SegNum,1)+CSI(SegNum,2)-1;
if numTRsVS ~= mr_ts_trs
	disp('TRs from Valid Segments txt and mgh file do not match. Fix it.')
	return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% uptake surface data
SubjectsFolder = '/oak/stanford/groups/leanew1/users/apines/surf';
% for surface data
surfL = [SubjectsFolder '/lh.sphere'];
surfR = [SubjectsFolder '/rh.sphere'];
% surface topography
[vx_l, faces_l] = read_surf(surfL);
[vx_r, faces_r] = read_surf(surfR);
% +1 the faces: begins indexing at 0
faces_l = faces_l + 1;
faces_r = faces_r + 1;
% normalize verts to unit sphere
% left
numV=length(vx_l);
vx_l(numV+1:end, :) = VecNormalize(vx_l(numV+1:end, :));
% right
numV=length(vx_r);
vx_r(numV+1:end, :) = VecNormalize(vx_r(numV+1:end, :));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% uptake functional data (on surface)
% handle input data
disp('Size of input data')
sizeInDl=size(TRs_l_g)
sizeInDr=size(TRs_r_g)
disp('Number of TRs detected L and R')
TR_n=sizeInDl(2)
sizeInDr(2)
assert(sizeInDl(2) == sizeInDr(2), 'Unequal time series length between hemispheres')

% left hemi
disp('converting left hemi to struct')
fl=struct;
% populate struct
for TRP=1:TR_n;
	fl.TRs{TRP}=TRs_l_g(:,TRP);
end

% r h
disp('converting right hemi to struct')
fr=struct;
for TRP=1:TR_n;
	fr.TRs{TRP}=TRs_r_g(:,TRP);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% compute optical flow on every pair of sequential TRs
% initialize output struct
us=struct;
disp('Computing optical flow: resting-state');

% initialize TRP counter: for plopping u outputs into master struct w/o/r/t their segment
% note trp = tr pair
TRPC=1;

% for each continuous segment
for seg=1:SegNum;
	% just to print out count of current segment
	seg
	SegStart=CSI(seg,1);
	SegSpan=CSI(seg,2);
	% get corresponding TRs from aggregate time series
	segTS_l=fl.TRs(SegStart:(SegStart+SegSpan-1));
	segTS_r=fr.TRs(SegStart:(SegStart+SegSpan-1));
	% loop over each TR-Pair: 1 fewer pair than number of TRs
	for TRP=1:(SegSpan-1)
		% print TR pair iter
		TRP
		% Compute decomposition.
		tic;
		% pull out adjacent frames
		u = of(N, faces_l, vx_l, segTS_l{TRP}, segTS_l{TRP+1}, h, alpha, s);
		% throw u into struct
		us.vf_left{TRPC}=u;
		% now right hemi
		u = of(N, faces_r, vx_r, segTS_r{TRP}, segTS_r{TRP+1}, h, alpha, s);
		toc;
		% throw u into struct
		us.vf_right{TRPC}=u;
		% update TR pair counter, which should increase +1 across segments
		TRPC=TRPC+1;
	end
end

save([childfp '/' subj '_' sesh '_' task '_OpFl.mat'],'us')


