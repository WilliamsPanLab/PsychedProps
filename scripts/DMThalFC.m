function ThalFC(subj)
% load in concatenated cifti, calculate FC for this subject, print out onto two maps (L and R)
Paths{1} = '/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti/tool_folder'; 
addpath(genpath(Paths{1}))
% set cifti path
ciftipath = ['/scratch/groups/leanew1/xcpd_outP50_36p_bp/xcp_d/' subj  '/ses-00/func/concatenated.dtseries.nii'];
% load in cifti
concatData=read_cifti(ciftipath);

% load in atlas
scale3fp='/oak/stanford/groups/leanew1/users/apines/maps/Tian_Subcortex_S3_3T_32k.dlabel.nii';
scale_3=read_cifti(scale3fp);
% load in scale 3 data
% 9 for right DM thal
DMT_rInd=find(scale_3.cdata==9);
% 34 for left DM thal
DMT_lInd=find(scale_3.cdata==34);

% extract DM thalamic TS - left and right
% left
DMT_ts_L=concatData.cdata(DMT_lInd,:);
DMT_mts_L=mean(DMT_ts_L);
% right
DMT_ts_R=concatData.cdata(DMT_rInd,:);
DMT_mts_R=mean(DMT_ts_R);

% get num of grayords
numGrayords=size(concatData.cdata);
numGrayords=numGrayords(1);
% initialize L amyg FC vec
LdmtFC=zeros(1,numGrayords);
% R
RdmtFC=zeros(1,numGrayords);
% get correlations across rest of brain (grayordinate-wise)
for g=1:numGrayords
	LdmtFC(g)=corr(concatData.cdata(g,:)',DMT_mts_L');
	RdmtFC(g)=corr(concatData.cdata(g,:)',DMT_mts_R');
end
% load in template cifti
tCifti = read_cifti('/oak/stanford/groups/leanew1/users/apines/maps/hcp.gradients.dscalar.nii')
% replace
tCifti.cdata(:,1)=LdmtFC;
tCifti.cdata(:,2)=RdmtFC;
% saveout into new cifti
write_cifti(tCifti,['/scratch/groups/leanew1/xcpd_outP50_36p_bp/bv_concat/' subj  '_dmtFC.dscalar.nii']);
