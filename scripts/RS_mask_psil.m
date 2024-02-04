function RS_mask_psil(subj,sesh)
% add paths
ToolFolder='/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti/tool_folder';
addpath(genpath(ToolFolder));

% load in boolean scan boundaries tmask to glean RS
Boundmaskfp=['/scratch/users/apines/PsiloData/' subj '/' subj '_' sesh '/func/run_boundaries_tmask.txt'];
Boundmask=load(Boundmaskfp);

% load in aggregate cifti 
CiftiFp=['/scratch/users/apines/PsiloData/' subj '/' subj '_' sesh '/func/' subj '_' sesh '_rsfMRI_uout_bpss_sr_noGSR_sm4.dtseries.nii'];
Cifti_file=read_cifti(CiftiFp);

% ensure length of boolean RS mask matches length of aggregate cifti
if size(Cifti_file.cdata, 2) ~= length(Boundmask)
    error('Error: Dimensions of boolean RS mask do not match the length of the aggregate CIFTI file.');
end

%%% extract resting-state portion of boundmask
% get run boundaries
runBounds=find(Boundmask==0);
% get run indices
runInds=bwlabel(Boundmask);
% get run 1
run1Ind=find(runInds==1);
% get run 2
run2Ind=find(runInds==2);
% ensure they are expected length
if length(run1Ind) ~=512
	if length(run1Ind) ~=509
		error(['unexpected RS length (rs1): ' num2str(length(run1Ind))])
	end
end
% mask to make only resting state
Cifti_rs1=Cifti_file.cdata(:,run1Ind);
% finagle cifti header to let me save the file
Cifti_file.diminfo{2}.length=length(run1Ind);
% implant rs1
Cifti_file.cdata=Cifti_rs1;
% save out cifti
CiftiFp=['/scratch/users/apines/PsiloData/' subj '/' subj '_' sesh '/func/' subj '_' sesh '_rs1.dtseries.nii']
write_cifti(Cifti_file,CiftiFp);

% for rs 2
if length(run2Ind) ~=512
	if length(run2Ind) ~=509
		error (['unexpected RS length (rs2): ' num2str(length(run2Ind))])
	end
end
% re-initialize cifti
CiftiFp=['/scratch/users/apines/PsiloData/' subj '/' subj '_' sesh '/func/' subj '_' sesh '_rsfMRI_uout_bpss_sr_noGSR_sm4.dtseries.nii'];
Cifti_file=read_cifti(CiftiFp);
% mask to make only resting state
Cifti_rs2=Cifti_file.cdata(:,run2Ind);
% finagle cifti header to let me save the file
Cifti_file.diminfo{2}.length=length(run2Ind);
% implant rs1
Cifti_file.cdata=Cifti_rs2;
% save out cifti
CiftiFp=['/scratch/users/apines/PsiloData/' subj '/' subj '_' sesh '/func/' subj '_' sesh '_rs2.dtseries.nii']
write_cifti(Cifti_file,CiftiFp);
