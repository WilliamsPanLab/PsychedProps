function Get_DMN_FC_human(subj) 
% add paths
addpath(genpath('/oak/stanford/groups/leanew1/users/apines/libs/'))
% ts filepath
% load in DMN
Dnet_32k=read_cifti('~/DMN_92k.dscalar.nii')';
% use DMN mask and fslr time series to get average DMN time series
childfp=['/scratch/groups/leanew1/xcpd_outP50_36p_bp/bv_rs/'];
fp=[childfp subj '_concatenated.dtseries.nii'];
data=read_cifti(fp);
TRs=data.cdata;
% make DMN mas
DMNbool=zeros(91282,1);
% 1:59412 to keep it cortical
DMNbool=Dnet_32k.cdata(1:59412,1)>.3;
% get avg DMN time series across both hemis
DMNTS=TRs(DMNbool,:);
meanDMNTS=mean(DMNTS);
% calculate average correlation with DMN time series
DMNcor = zeros(size(TRs, 1), 1);
% for each grayordinate
for i = 1:size(TRs, 1)
    DMNcor(i) = corr(meanDMNTS', TRs(i, :)');
end
%output cifti
template_cifti = Dnet_32k;
template_cifti.cdata = DMNcor;
outFP=[childfp '/' subj '_DMNFC_map.dscalar.nii'];
ciftisave(template_cifti,outFP);
% make output structure
CorStruct={};
CorStruct.DMNcor=DMNcor;
% save out as a 1xn vertices vector for each subject with subj sesh task in fn for merging later
save([childfp '/' subj '_DMNFC_map.mat'],'CorStruct');
