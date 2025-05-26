function RS_mask_psil(subj,sesh)
% add paths
ToolFolder='/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti/tool_folder';
addpath(genpath(ToolFolder));

% load in boolean scan boundaries tmask to glean RS
Boundmaskfp=['/oak/stanford/groups/leanew1/SHARED_DATASETS/private/WashU_psilocybin/' subj '/' subj '_' sesh '/func/run_boundaries_tmask.txt'];
Boundmask=load(Boundmaskfp);

% load in aggregate cifti 
CiftiFp=['/oak/stanford/groups/leanew1/SHARED_DATASETS/private/WashU_psilocybin/' subj '/' subj '_' sesh '/func/' subj '_' sesh '_rsfMRI_uout_bpss_sr_noGSR_sm4.dtseries.nii'];
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

% get number of RS scans
tabulatedRunsInds=tabulate(runInds);
LongScans=find(tabulatedRunsInds(:,2)>500);
LongScanInds=tabulatedRunsInds(LongScans,1);
% print number of rs scans found
disp([num2str(length(LongScanInds)) ' rs scans found for ' subj ' session ' sesh])

% alter to run for each RS scan
for run=1:length(LongScanInds)
	% re-initialize cifti
	CiftiFp=['/oak/stanford/groups/leanew1/SHARED_DATASETS/private/WashU_psilocybin/' subj '/' subj '_' sesh '/func/' subj '_' sesh '_rsfMRI_uout_bpss_sr_noGSR_sm4.dtseries.nii'];
	Cifti_file=read_cifti(CiftiFp);
	runInd=find(runInds==LongScanInds(run));
	% ensure expected length
	if length(runInd) ~=512
   	     	if length(runInd) ~=509
               		error(['unexpected RS length run' num2str(LongScanInds(run)) ': ' num2str(length(runInd))])
        	end
	end
	% mask master cifti for resting-only saveout
	Cifti_rs=Cifti_file.cdata(:,runInd);
	% implant masked time series
	Cifti_file.cdata=Cifti_rs;
	Cifti_file.diminfo{2}.length=length(runInd);
	CiftiFp=['/oak/stanford/groups/leanew1/SHARED_DATASETS/private/WashU_psilocybin/' subj '/' subj '_' sesh '/func/' subj '_' sesh '_rs' num2str(run) '.dtseries.nii']
	write_cifti(Cifti_file,CiftiFp);
end	
