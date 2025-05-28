function SubCort_OpFl_caud_R(subj,sesh,task)
% all credit to NeuroPattToolbox: https://github.com/BrainDynamicsUSYD/NeuroPattToolbox, at least Rory Townsend, Xian Long, and Pulin Gong
addpath(genpath('/oak/stanford/groups/leanew1/users/apines/scripts/NeuroPattToolbox'));
addpath(genpath('/oak/stanford/groups/leanew1/users/apines/libs/'))
% Most of code is from:
% Rory Townsend, Aug 2018
% rory.townsend@sydney.edu.au
% set sampling frequency
Fs=(1/.71);
params=struct;
params.zscoreChannels=0;
params.params.subtractBaseline=0;
% enforcement of smoothness of vector field
params.opAlpha = 0.5;
startTime=datetime
% https://github.com/rorygt/NeuroPattToolbox/blob/master/setParams.m claims opBeta=0.01 is default, but it's not clear we can run the code with that low of a beta and 1 is what is actually set in the code. Trying 0.1 to afford some nonlinearity while being computationally tractable
% enforcement of linearity of vector field (higher = more linear)
params.opBeta = 0.1;
% note we don't want to analyze phase directly in case propagations of interest are aperiodic
params.useAmplitude = true;

% load in time series for this subj/sesh
childfp=['/scratch/groups/leanew1/xcpd_outP50_36p_bp/xcp_d/' subj '/' sesh '/func' ];

% load in data
if task=='rs1'
	fp=[childfp '/' subj '_' sesh '_task-rs_acq-mb_dir-pe0_run-0_space-fsLR_den-91k_desc-denoisedSmoothed_bold.dtseries.nii'];
elseif task=='rs2'
	fp=[childfp '/' subj '_' sesh '_task-rs_acq-mb_dir-pe1_run-0_space-fsLR_den-91k_desc-denoisedSmoothed_bold.dtseries.nii'];
else
	fp=[childfp '/' subj '_' sesh '_task-' task '_acq-mb_dir-pe0_run-0_space-fsLR_den-91k_desc-denoisedSmoothed_bold.dtseries.nii'];
end
% convert it to a nifti
niftiout = [childfp '/' subj '_' sesh '_task_' task '_volumetric_subcort.nii.gz'];
command = ['wb_command -cifti-separate ' fp ' COLUMN -volume-all ' niftiout];
system(command)
% to let system catch up
pause(5);
data=niftiread(niftiout);
% read in temporal mask
subjDir=['/scratch/users/apines/data/mdma/' subj '/' sesh];
tempMask=[subjDir '/' subj '_' sesh '_task-' task '_AllSegments.txt'];
tempMask=dlmread(tempMask);
% read in atlas
atlas=read_cifti('/oak/stanford/groups/leanew1/users/apines/maps/Tian_Subcortex_S1_3T_32k.dscalar.nii');
% extract voxels corresponding to subcortex
subcortVoxels=zeros(0,3);
% first 2 are cortex
for s=3:length(atlas.diminfo{1}.models);
	voxList=atlas.diminfo{1}.models{s}.voxlist;
	subcortVoxels=vertcat(subcortVoxels,voxList);
end
% now get where in this list subcort structures of interst are
atlasInds=atlas.cdata;
%https://github.com/yetianmed/subcortex/blob/62751401be0dba587442799748a886d5d87b8364/Group-Parcellation/3T/Subcortex-Only/Tian_Subcortex_S1_3T_label.txt
% segment out caudate: numeric labels are 8 for right and 16 for left
% need to subtract off number of cortical verts to match subcortVoxel list
caud_R=find(atlasInds==8)-59412;
caud_L=find(atlasInds==16)-59412;
% segment out hippocampus: numeric labels are 1 for right and 9 for left
hippo_R=find(atlasInds==1)-59412;
hippo_L=find(atlasInds==9)-59412;
% pull actual coordinates: zero-indexing correction with +1 
caud_R_coords=subcortVoxels(caud_R,:)+1;
caud_L_coords=subcortVoxels(caud_L,:)+1;
hippo_R_coords=subcortVoxels(hippo_R,:)+1;
hippo_L_coords=subcortVoxels(hippo_L,:)+1;
% extract boxes that these voxels live in (we're going to restrict box by 1 to live less on the edge)
caud_R_min = min(caud_R_coords,[],1)+1;
caud_R_max = max(caud_R_coords,[],1)-1;
% caud L
caud_L_min = min(caud_L_coords,[],1)+1;
caud_L_max = max(caud_L_coords,[],1)-1;
% hippo R
hippo_R_min = min(hippo_R_coords,[],1)+1;
hippo_R_max = max(hippo_R_coords,[],1)-1;
% left
hippo_L_min = min(hippo_L_coords,[],1)+1;
hippo_L_max = max(hippo_L_coords,[],1)-1;
% caud R min needs an additional 1-pixel shave-off: I will attempt to depict why visually.
% Box outlining Left caudate min/max |--/----|
%                                    | /     |
%          caudate as "/"            |/------|
%                                    /
% in words, the -1 to boundaries we apply to minimize boundary effects on optical flow in the z direction reduces the real y-extent of caudate signal. This is because of the diagonal entry of the caudate into our bounding box.
% so with conservative z-boundaries, the lowest-numbered y-plane is actually all 0
caud_R_min(2)=caud_R_min(2)+1;

% set caud R mask
caud_R_mask = false(size(data,1), size(data,2), size(data,3));
% now set voxels in caud R as true
for i = 1:size(caud_R_coords,1)
    x = caud_R_coords(i,1);
    y = caud_R_coords(i,2);
    z = caud_R_coords(i,3);
    caud_R_mask(x, y, z) = true;
end
% now apply this mask to data
for t = 1:size(data,4)
    vol = data(:,:,:,t);
    vol(~caud_R_mask) = 0;
    data(:,:,:,t) = vol;
end

% now zoom in on the bounding box
subcortVol = data(caud_R_min(1):caud_R_max(1), ...
                caud_R_min(2):caud_R_max(2), ...
                caud_R_min(3):caud_R_max(3), :);

% visualize some timepoints if of interest
%timepoints = 72:82;
% arbitrary slice indices
%sliInd=3;
%outdir = './hippo_L_masked_bounded_pngs/';
%if ~exist(outdir, 'dir'), mkdir(outdir); end
%for t = timepoints
%    vol = double(data(hippo_L_min(1):hippo_L_max(1), ...
%                      hippo_L_min(2):hippo_L_max(2), ...
%                      hippo_L_min(3):hippo_L_max(3), t));
    % mask out what is not left hippocampus
%    vol(~hippo_L_mask) = NaN;
%    slice = vol(:,:,sliInd);  % axial slice at z = sliInd
    % this accounts for difference in orientation and what matlab depicts
%    oriented_slice = rot90(flipud(slice));  % flip Y, rotate CCW
%    f = figure('visible', 'off');
%    imagesc(oriented_slice, [0 20]);
%    colormap(jet); axis image off; colorbar;
%    title(sprintf('hippo_L | t=%d | z=%d', t, sliInd));
%    print(f, fullfile(outdir, sprintf('hippoL_t%03d_z%02d.png', t, sliInd)), '-dpng', '-r300');
%    close(f);
%end

% apply temporal mask to subcortVol
tempMaskIncl=tempMask((tempMask(:,3)==1),:);
% count number of good segments
goodSegs=size(tempMaskIncl,1);
% get boolean indices of timepoints to include
booleanInds=false(1,size(subcortVol,4));
for i = 1:goodSegs
        booleanInds(tempMaskIncl(i,1):tempMaskIncl(i,2)) = true;
end
% mask the time series
subcortVol=subcortVol(:,:,:,booleanInds);

% run on each slice in Z
[nx, ny, nz, nt] = size(subcortVol);
% pre-allocate vfs
vfs_x = zeros(size(subcortVol));
vfs_x = vfs_x(:,:,:,1:end-1);
vfs_y = zeros(size(subcortVol));
vfs_y = vfs_y(:,:,:,1:end-1);
vfs_z = zeros(size(subcortVol));
vfs_z = vfs_z(:,:,:,1:end-1);
% mask
mask_cropped = all(subcortVol ~= 0, 4);
% maintaining their variable names
badChannels3D = logical(~mask_cropped);
% initialize tr pair counter
trpc=1;
% for each segment
for s=1:goodSegs
s
SegStart=tempMaskIncl(s,1);
SegSpan=tempMaskIncl(s,4);
segTS=subcortVol(:,:,:,trpc:(trpc+SegSpan-1));
% for each x slice
for x = 1:nx
        wvcfs = subcortVol(x,:,:,:);
        badChannels = badChannels3D(x,:,:);
	% squeeze so we are in dim 1 x dim 2 x time
        wvcfs2D = double(squeeze(wvcfs));
	badChannels=squeeze(badChannels);
        % run optical flow
        [vx, vy, csteps] = opticalFlow2(wvcfs2D, badChannels, ...
            params.opAlpha, params.opBeta, ~params.useAmplitude);
	% populate it: -2 to combine with -1 offset from above (inclusive indexing), -1 more because opflow is on tr pairs not trs
         vfs_x(x,:,:,trpc:(trpc+SegSpan-2)) = vx + 1i*vy;
end
% for each y slice
for y = 1:ny
        wvcfs = subcortVol(:,y,:,:);
        badChannels = badChannels3D(:,y,:);
        % squeeze so we are in dim 1 x dim 2 x time
        wvcfs2D = double(squeeze(wvcfs));
	badChannels=squeeze(badChannels);
        % run optical flow
        [vx, vy, csteps] = opticalFlow2(wvcfs2D, badChannels, ...
            params.opAlpha, params.opBeta, ~params.useAmplitude);
	% populate it: -2 to combine with -1 offset from above (inclusive indexing), -1 more because opflow is on tr pairs not trs
         vfs_y(:,y,:,trpc:(trpc+SegSpan-2)) = vx + 1i*vy;
end
% for each z slice
for z = 1:nz
	wvcfs = subcortVol(:,:,z,:);
        badChannels = badChannels3D(:,:,z);
	% squeeze so we are in dim 1 x dim 2 x time
	wvcfs2D = double(squeeze(wvcfs));
	badChannels=squeeze(badChannels);
	% run optical flow
	[vx, vy, csteps] = opticalFlow2(wvcfs2D, badChannels, ...
            params.opAlpha, params.opBeta, ~params.useAmplitude);
	% populate it: -2 to combine with -1 offset from above (inclusive indexing), -1 more because opflow is on tr pairs not trs
         vfs_z(:,:,z,trpc:(trpc+SegSpan-2)) = vx + 1i*vy;
end
trpc=trpc+SegSpan;
end

% Set all outputs
% AP omitting onlyPatterns: we want dem vecta fieldz
outputs.velocityFields_x = vfs_x;
outputs.velocityFields_y = vfs_y;
outputs.velocityFields_z = vfs_z;
outputs.badChannels = badChannels;
outputs.Fs = Fs;
outputs.processTime = datetime - startTime;

% save outputs
outfp=[childfp '/' subj '_CaudR_vf_out_' sesh '_' task '.mat'];
save(outfp,'outputs');
