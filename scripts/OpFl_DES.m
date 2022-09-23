function us = OpFl_DES(subj)
% add paths
addpath(genpath('/oak/stanford/groups/leanew1/users/apines/libs/'))
%%%%%%%%%%%%%%%%%%%% Set parameters
N = 10; % Vector spherical harmonics degree
h = 1; % finite-difference height vector of triangles
alpha = 1; % scales Tikhonov regularization, > alpha = > spatiotemporal smoothness
s = 1; % R(u), regularizing functional, scales Tikhonov regularization more rapidly via penalization of first spatial and first temporal derivative
%%%%%%%%%%%%%%%%%%%%

% load in data
fpL=['/scratch/users/apines/data/DES/' subj '/wm_lh_cent.mgh']
fpR=['/scratch/users/apines/data/DES/' subj '/wm_rh_cent.mgh']
dataL=MRIread(fpL).vol;
dataR=MRIread(fpR).vol;
% squeeze to get rid of extra dimensions
TRs_l_w=squeeze(dataL);
TRs_r_w=squeeze(dataR);
% noncon
fpL=['/scratch/users/apines/data/DES/' subj '/rs_lh_cent.mgh']
fpR=['/scratch/users/apines/data/DES/' subj '/rs_rh_cent.mgh']
dataL=MRIread(fpL).vol;
dataR=MRIread(fpR).vol;
TRs_l_rs=squeeze(dataL);
TRs_r_rs=squeeze(dataR);
% con
fpL=['/scratch/users/apines/data/DES/' subj '/rs2_lh_cent.mgh']
fpR=['/scratch/users/apines/data/DES/' subj '/rs2_rh_cent.mgh']
dataL=MRIread(fpL).vol;
dataR=MRIread(fpR).vol;
TRs_l_rs2=squeeze(dataL);
TRs_r_rs2=squeeze(dataR);

% LEGACY COMMENTED OUT UNTIL DATA PROC'ED
% fpL=['/scratch/users/apines/data/ispot/' subj '/gng_lh_cent.mgh']
% fpR=['/scratch/users/apines/data/ispot/' subj '/gng_rh_cent.mgh']
% dataL=MRIread(fpL).vol;
% dataR=MRIread(fpR).vol;
% squeeze to get rid of extra dimensions
% TRs_l_g=squeeze(dataL);
% TRs_r_g=squeeze(dataR);
% noncon
% fpL=['/scratch/users/apines/data/ispot/' subj '/noncon_lh_cent.mgh']
% fpR=['/scratch/users/apines/data/ispot/' subj '/noncon_rh_cent.mgh']
% dataL=MRIread(fpL).vol;
% dataR=MRIread(fpR).vol;
% TRs_l_ncf=squeeze(dataL);
% TRs_r_ncf=squeeze(dataR);
% con
% fpL=['/scratch/users/apines/data/ispot/' subj '/con_lh_cent.mgh']
% fpR=['/scratch/users/apines/data/ispot/' subj '/con_rh_cent.mgh']
% dataL=MRIread(fpL).vol;
% dataR=MRIread(fpR).vol;
% TRs_l_con=squeeze(dataL);
% TRs_r_con=squeeze(dataR);

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% WM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% uptake functional data (on surface)
% handle input data
disp('Size of input data')
sizeInDl=size(TRs_l_w)
sizeInDr=size(TRs_r_w)
disp('Number of TRs detected L and R')
TR_n=sizeInDl(2)
sizeInDr(2)
assert(sizeInDl(2) == sizeInDr(2), 'Unequal time series length between hemispheres')

% left hemi
disp('converting left hemi to struct')
fl=struct;
% populate struct
for TRP=1:TR_n;
        fl.TRs{TRP}=TRs_l_w(:,TRP);
end

% r h
disp('converting right hemi to struct')
fr=struct;
for TRP=1:TR_n;
        fr.TRs{TRP}=TRs_r_w(:,TRP);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% compute optical flow on every pair of sequential TRs
% initialize output struct
us=struct;
disp('Computing optical flow: Continuous Performance');

% initialize TRP counter: for plopping u outputs into master struct w/o/r/t their segment
% note trp = tr pair
TRPC=1;

% loop over each TR-Pair: 1 fewer pair than number of TRs
for TRP=1:(TR_n-1)
        % print TR pair iter
        TRP
        % Compute decomposition.
        tic;
        % pull out adjacent frames
        u = of(N, faces_l, vx_l, fl.TRs{TRP}, fl.TRs{TRP+1}, h, alpha, s);
        % throw u into struct
        us.vf_left{TRPC}=u;
        % now right hemi
        u = of(N, faces_r, vx_r, fr.TRs{TRP}, fr.TRs{TRP+1}, h, alpha, s);
        toc;
        % throw u into struct
        us.vf_right{TRPC}=u;
        % update TR pair counter, which should increase +1 across segments
        TRPC=TRPC+1;
end

save(['/scratch/users/apines/data/DES/' subj '/OpFl_wm.mat'],'us')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RS 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% uptake functional data (on surface)
% handle input data
disp('Size of input data')
sizeInDl=size(TRs_l_rs)
sizeInDr=size(TRs_r_rs)
disp('Number of TRs detected L and R')
TR_n=sizeInDl(2)
sizeInDr(2)
assert(sizeInDl(2) == sizeInDr(2), 'Unequal time series length between hemispheres')

% left hemi
disp('converting left hemi to struct')
fl=struct;
% populate struct
for TRP=1:TR_n;
        fl.TRs{TRP}=TRs_l_rs(:,TRP);
end

% r h
disp('converting right hemi to struct')
fr=struct;
for TRP=1:TR_n;
        fr.TRs{TRP}=TRs_r_rs(:,TRP);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% compute optical flow on every pair of sequential TRs
% initialize output struct
us=struct;
disp('Computing optical flow: Resting State');

% initialize TRP counter: for plopping u outputs into master struct w/o/r/t their segment
% note trp = tr pair
TRPC=1;

% loop over each TR-Pair: 1 fewer pair than number of TRs
for TRP=1:(TR_n-1)
        % print TR pair iter
        TRP
        % Compute decomposition.
        tic;
        % pull out adjacent frames
        u = of(N, faces_l, vx_l, fl.TRs{TRP}, fl.TRs{TRP+1}, h, alpha, s);
        % throw u into struct
        us.vf_left{TRPC}=u;
        % now right hemi
        u = of(N, faces_r, vx_r, fr.TRs{TRP}, fr.TRs{TRP+1}, h, alpha, s);
        toc;
        % throw u into struct
        us.vf_right{TRPC}=u;
        % update TR pair counter, which should increase +1 across segments
        TRPC=TRPC+1;
end

save(['/scratch/users/apines/data/ispot/' subj '/OpFl_rs.mat'],'us')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RS 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% uptake functional data (on surface)
% handle input data
disp('Size of input data')
sizeInDl=size(TRs_l_rs2)
sizeInDr=size(TRs_r_rs2)
disp('Number of TRs detected L and R')
TR_n=sizeInDl(2)
sizeInDr(2)
assert(sizeInDl(2) == sizeInDr(2), 'Unequal time series length between hemispheres')

% left hemi
disp('converting left hemi to struct')
fl=struct;
% populate struct
for TRP=1:TR_n;
        fl.TRs{TRP}=TRs_l_rs2(:,TRP);
end

% r h
disp('converting right hemi to struct')
fr=struct;
for TRP=1:TR_n;
        fr.TRs{TRP}=TRs_r_rs2(:,TRP);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% compute optical flow on every pair of sequential TRs
% initialize output struct
us=struct;
disp('Computing optical flow: Resting State 2');

% initialize TRP counter: for plopping u outputs into master struct w/o/r/t their segment
% note trp = tr pair
TRPC=1;

% loop over each TR-Pair: 1 fewer pair than number of TRs
for TRP=1:(TR_n-1)
        % print TR pair iter
        TRP
        % Compute decomposition.
        tic;
        % pull out adjacent frames
        u = of(N, faces_l, vx_l, fl.TRs{TRP}, fl.TRs{TRP+1}, h, alpha, s);
        % throw u into struct
        us.vf_left{TRPC}=u;
        % now right hemi
        u = of(N, faces_r, vx_r, fr.TRs{TRP}, fr.TRs{TRP+1}, h, alpha, s);
        toc;
        % throw u into struct
        us.vf_right{TRPC}=u;
        % update TR pair counter, which should increase +1 across segments
        TRPC=TRPC+1;
end

save(['/scratch/users/apines/data/ispot/' subj '/OpFl_rs2.mat'],'us')

%
%%%%%
%%%%%%%%%%
%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LEGACY
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%
%%%%%%%%%%
%%%%%
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GO NO GO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
disp('Computing optical flow: Go No Go');

% initialize TRP counter: for plopping u outputs into master struct w/o/r/t their segment
% note trp = tr pair
TRPC=1;

% loop over each TR-Pair: 1 fewer pair than number of TRs
for TRP=1:(TR_n-1)
	% print TR pair iter
	TRP
	% Compute decomposition.
	tic;
	% pull out adjacent frames
	u = of(N, faces_l, vx_l, fl.TRs{TRP}, fl.TRs{TRP+1}, h, alpha, s);
	% throw u into struct
	us.vf_left{TRPC}=u;
	% now right hemi
	u = of(N, faces_r, vx_r, fr.TRs{TRP}, fr.TRs{TRP+1}, h, alpha, s);
	toc;
	% throw u into struct
	us.vf_right{TRPC}=u;
	% update TR pair counter, which should increase +1 across segments
	TRPC=TRPC+1;
end

save(['/scratch/users/apines/data/ispot/' subj '/OpFl_GNG.mat'],'us')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NONCON FACES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% uptake functional data (on surface)
% handle input data
disp('Size of input data')
sizeInDl=size(TRs_l_ncf)
sizeInDr=size(TRs_r_ncf)
disp('Number of TRs detected L and R')
TR_n=sizeInDl(2)
sizeInDr(2)
assert(sizeInDl(2) == sizeInDr(2), 'Unequal time series length between hemispheres')

% left hemi
disp('converting left hemi to struct')
fl=struct;
% populate struct
for TRP=1:TR_n;
        fl.TRs{TRP}=TRs_l_ncf(:,TRP);
end

% r h
disp('converting right hemi to struct')
fr=struct;
for TRP=1:TR_n;
        fr.TRs{TRP}=TRs_r_ncf(:,TRP);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% compute optical flow on every pair of sequential TRs
% initialize output struct
us=struct;
disp('Computing optical flow: NonCon Faces');

% initialize TRP counter: for plopping u outputs into master struct w/o/r/t their segment
% note trp = tr pair
TRPC=1;

% loop over each TR-Pair: 1 fewer pair than number of TRs
for TRP=1:(TR_n-1)
        % print TR pair iter
        TRP
        % Compute decomposition.
        tic;
        % pull out adjacent frames
        u = of(N, faces_l, vx_l, fl.TRs{TRP}, fl.TRs{TRP+1}, h, alpha, s);
        % throw u into struct
        us.vf_left{TRPC}=u;
        % now right hemi
        u = of(N, faces_r, vx_r, fr.TRs{TRP}, fr.TRs{TRP+1}, h, alpha, s);
        toc;
        % throw u into struct
        us.vf_right{TRPC}=u;
        % update TR pair counter, which should increase +1 across segments
        TRPC=TRPC+1;
end

save(['/scratch/users/apines/data/ispot/' subj '/OpFl_NCF.mat'],'us')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CON FACES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% uptake functional data (on surface)
% handle input data
disp('Size of input data')
sizeInDl=size(TRs_l_con)
sizeInDr=size(TRs_r_con)
disp('Number of TRs detected L and R')
TR_n=sizeInDl(2)
sizeInDr(2)
assert(sizeInDl(2) == sizeInDr(2), 'Unequal time series length between hemispheres')

% left hemi
disp('converting left hemi to struct')
fl=struct;
% populate struct
for TRP=1:TR_n;
        fl.TRs{TRP}=TRs_l_con(:,TRP);
end

% r h
disp('converting right hemi to struct')
fr=struct;
for TRP=1:TR_n;
        fr.TRs{TRP}=TRs_r_con(:,TRP);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% compute optical flow on every pair of sequential TRs
% initialize output struct
us=struct;
disp('Computing optical flow: Con Faces');

% initialize TRP counter: for plopping u outputs into master struct w/o/r/t their segment
% note trp = tr pair
TRPC=1;

% loop over each TR-Pair: 1 fewer pair than number of TRs
for TRP=1:(TR_n-1)
        % print TR pair iter
        TRP
        % Compute decomposition.
        tic;
        % pull out adjacent frames
        u = of(N, faces_l, vx_l, fl.TRs{TRP}, fl.TRs{TRP+1}, h, alpha, s);
        % throw u into struct
        us.vf_left{TRPC}=u;
        % now right hemi
        u = of(N, faces_r, vx_r, fr.TRs{TRP}, fr.TRs{TRP+1}, h, alpha, s);
        toc;
        % throw u into struct
        us.vf_right{TRPC}=u;
        % update TR pair counter, which should increase +1 across segments
        TRPC=TRPC+1;
end

save(['/scratch/users/apines/data/ispot/' subj '/OpFl_CON.mat'],'us')
