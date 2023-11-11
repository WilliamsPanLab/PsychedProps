function Extract_NGSC(subj,sesh,task)
% set parent directory
parentfp=['/scratch/groups/leanew1/xcpd_outP50_36p_bp/xcp_d/' subj '/' sesh '/func'];

% define some paths 
Paths{1} = '/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti/tool_folder';
addpath(genpath(Paths{1}))
% read in time series
C=ft_read_cifti_mod([parentfp '/' subj '_' sesh '_task-' task '_acq-mb_dir-pe0_run-0_space-fsLR_den-91k_desc-denoised_bold.dtseries.nii']; 
timeseries=C.data;

% load in default mode network boundaries
dmn=ft_read_cifti_mod('/oak/stanford/groups/leanew1/users/apines/data/RobustInitialization/SoftParcel.dscalar.nii');
% mask as 0.3 or greater
dmn.data(dmn.data<0.3)=0;
dmn.data(dmn.data>=0.3)=1;
% apply to time series
dmn_ts=timeseries(dmn.data==1,:);

% conduct temporal PCA to recover components m (# of grayordinates)
[coeff,score,latent,tsquared,explained,mu] = pca(dmn_ts);

% get normalized eigenvalue of each principal component (divided by sum of eigenvalues)
norm_eig=explained/sum(explained);?

% get the numerator for normalized entropy of the normalized eigenvalues, sum of the normalized eigenvalue * log(normalized eigenvalue) for each component

% get the denominator for normalized entropy of the normalized eigenvalues, log of the number of components

% get normalized entropy of the normalized eigenvalues, numerator divided by denominator * -1

% save out normalized entropy for dmn
