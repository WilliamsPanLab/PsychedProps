function Get_Thr_FC_human(subj) 
% add paths
addpath(genpath('/oak/stanford/groups/leanew1/users/apines/libs/'))
% ts filepath
childfp=['/scratch/groups/leanew1/xcpd_outP50_36p_bp/bv_rs/'];
fp=[childfp subj '_concatenated.dtseries.nii'];
data=read_cifti(fp);
TRs=data.cdata;

% Load Brodmann labels
Bro_L = gifti('/scratch/users/apines/lh.Brodmann.32k_fs_LR.label.gii');
Bro_R = gifti('/scratch/users/apines/rh.Brodmann.32k_fs_LR.label.gii');

% Load dtseries
data = read_cifti(fp);
TRs = data.cdata;
models = data.diminfo{1}.models;

% Get valid surface vertex indices used in CIFTI
vl_L = models{1}.vertlist+1;        % 0-based indexing
vl_R = models{2}.vertlist+1;        % 0-based indexing
start_R = models{2}.start;        % CIFTI index (starts at 29697)

% Initialize mask
combinedMask = false(91282, 1);

% get BA25 in nonmw indices
NoMWL=Bro_L.cdata(vl_L);
NoMWR=Bro_R.cdata(vl_R);
% note brod 25 is acutally 21 -_-
BA25_L=NoMWL==21;
BA25_R=NoMWR==21;

% Apply BA25 cortical mask
combinedMask(1:29696) = BA25_L;
combinedMask(29697:59412) = BA25_R;

% Add subcortical amygdala voxels
amyL = models{5};  % AMYGDALA_LEFT
amyR = models{6};  % AMYGDALA_RIGHT
combinedMask(amyL.start : amyL.start + amyL.count - 1) = true;
combinedMask(amyR.start : amyR.start + amyR.count - 1) = true;

% get avg threat time series across both hemis
TTS=TRs(combinedMask,:);
meanTTS=mean(TTS);
% calculate average correlation with DMN time series
Tcor = zeros(size(TRs, 1), 1);
% for each grayordinate
for i = 1:size(TRs, 1)
    Tcor(i) = corr(meanTTS', TRs(i, :)');
end
%output cifti
Dnet_32k=read_cifti('~/DMN_92k.dscalar.nii')';
template_cifti = Dnet_32k;
template_cifti.cdata = Tcor;
outFP=[childfp '/' subj '_TFC_map.dscalar.nii'];
ciftisave(template_cifti,outFP);
% make output structure
CorStruct={};
CorStruct.Tcor=Tcor;
% save out as a 1xn vertices vector for each subject with subj sesh task in fn for merging later
save([childfp '/' subj '_TFC_map.mat'],'CorStruct');
