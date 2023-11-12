function Extract_NGSC(subj,sesh,task)
% set parent directory
parentfp=['/scratch/groups/leanew1/xcpd_outP50_36p_bp/xcp_d/' subj '/' sesh '/func'];

% define some paths 
Paths{1} = '/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti/tool_folder';
addpath(genpath(Paths{1}))

% read in time series
C=ft_read_cifti_mod([parentfp '/' subj '_' sesh '_task-' task '_acq-mb_dir-pe0_run-0_space-fsLR_den-91k_desc-denoisedSmoothed_bold.dtseries.nii']); 
C_timeseries=C.data;

% read in gordon parcellation
GP=ft_read_cifti_mod('~/null_lL_WG33/Gordon333_FreesurferSubcortical.32k_fs_LR.dlabel.nii');

% gordon parcel labels
% https://balsa.wustl.edu/file/JX5V

% load in default mode network boundaries
dmn=ft_read_cifti_mod('/oak/stanford/groups/leanew1/users/apines/data/RobustInitialization/SoftParcel.dscalar.nii');
% second is DMN
dmn.data=dmn.data(:,2);

% get average dmn value for each parcel
pDMN=zeros(1,333);
for p=1:333
	pDMN(p)=mean(dmn.data(GP.data==p));
end
% get gordon parcels where average DMN value is .3 or greater
DMNParcels=find(pDMN>.3);
cxDMN=[];
% for each parcel, get complexity of timeseries
for p=DMNParcels
	% apply to time series
	dmn_ts=C_timeseries(logical(GP.data==p),:);
	% conduct temporal PCA to recover components m (# of grayordinates)
	[coeff, ~, explained] = pca(dmn_ts);
	%[coeff,score,latent,tsquared,explained,mu] = pca(dmn_ts);
	% get normalized eigenvalue of each principal component (divided by sum of eigenvalues)
	% https://stackoverflow.com/questions/30792185/matlab-how-to-obtain-the-eigenvalues-from-the-pca
	norm_eig=explained/sum(explained);

	% get the numerator for normalized entropy of the normalized eigenvalues, sum of the normalized eigenvalue * log(normalized eigenvalue) for each component
	numerator = sum(norm_eig .* log(norm_eig));

	% get the denominator for normalized entropy of the normalized eigenvalues, log of the number of components
	denominator = log(length(norm_eig));

	% get normalized entropy of the normalized eigenvalues, numerator divided by denominator * -1
	nGSC = -1 * (numerator / denominator);
	cxDMN = [cxDMN nGSC];
end
% can insert parcelwise saveout here later if interested
FullParcels=zeros(1,333);
FullParcels(DMNParcels)=1;
GPdataOut=zeros(91282,1);
iterator=1
for p=DMNParcels
   	GPdataOut(GP.data==p)=cxDMN(iterator);	
	iterator=iterator+1;
end
% save out
ComplOut = read_cifti('/oak/stanford/groups/leanew1/users/apines/NeuroSynthMaps/dmn_smooth.dscalar.nii')
ComplOut.cdata=GPdataOut;
ComplOut.diminfo{2}.length=1
ComplOut.diminfo{2}.type='series'
ComplOut.diminfo{2}.seriesUnit='SECOND'
ComplOut.diminfo{2}.seriesStep=.71
ComplOut.diminfo{2}.seriesStart=0.00
write_cifti(ComplOut,[parentfp '/' subj '_' sesh '_' task '_Complexity.dtseries.nii']);

% save out normalized entropy for dmn
avComplexity=mean(cxDMN);

T=table(avComplexity,'RowNames',"Row1");
outFP=['/scratch/users/apines/data/mdma/' subj '/' sesh];
writetable(T,[outFP '/' subj '_' sesh '_' task '_Complexity_gro.csv'],'WriteRowNames',true)


