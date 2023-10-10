function AmygFC(subj,sesh)
% load in concatenated cifti, calculate FC for this subject, print out onto two maps (L and R)
Paths{1} = '/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti/tool_folder'; 
addpath(genpath(Paths{1}))
% set cifti path
ciftipath = ['/scratch/groups/leanew1/xcpd_outP50_36p_bp/bv_concat/' subj  '-concat.dtseries.nii'];
% load in cifti
concatData=read_cifti(ciftipath);
% get left and right amygdala grayords/voxels
LAmygStart=concatData.diminfo{1}.models{5}.start;
RAmygStart=concatData.diminfo{1}.models{6}.start;
LAmygCount=concatData.diminfo{1}.models{5}.count;
RAmygCount=concatData.diminfo{1}.models{6}.count;
% extract amygdalar TS - left and right
% left
A_ts_L=concatData.cdata(LAmygStart:(LAmygStart+LAmygCount-1),:);
A_mts_L=mean(A_ts_L);
% right
A_ts_R=concatData.cdata(RAmygStart:(RAmygStart+RAmygCount-1),:);
A_mts_R=mean(A_ts_R);
% get num of grayords
numGrayords=size(concatData.cdata);
numGrayords=numGrayords(1);
% initialize L amyg FC vec
LamygFC=zeros(1,numGrayords);
% R
RamygFC=zeros(1,numGrayords);
% get correlations across rest of brain (grayordinate-wise)
for g=1:numGrayords
	LamygFC(g)=corr(concatData.cdata(g,:)',A_mts_L');
	RamygFC(g)=corr(concatData.cdata(g,:)',A_mts_R');
end
% load in template cifti
tCifti = read_cifti('/oak/stanford/groups/leanew1/users/apines/maps/hcp.gradients.dscalar.nii')
% replace
tCifti.cdata(:,1)=LamygFC;
tCifti.cdata(:,2)=RamygFC;
% saveout into new cifti
write_cifti(tCifti,['/scratch/groups/leanew1/xcpd_outP50_36p_bp/bv_concat/' subj  '_amygFC.dscalar.nii']);

