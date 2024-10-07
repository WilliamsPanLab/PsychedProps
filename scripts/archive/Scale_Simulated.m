function Scale_Simulated(subj,sesh,task)
% scale simulated data to be same magnitude as real data (range of signal)
addpath(genpath('/oak/stanford/groups/leanew1/users/apines/libs/'))
% read in sample data
parentfp=['/scratch/groups/leanew1/xcpd_outP50_36p_bp/xcp_d/' subj '/' sesh '/func'];
if string(task)=="rs1"
        fp=[parentfp '/' subj '_' sesh '_task-rs_acq-mb_dir-pe0_run-0_space-fsLR_den-91k_desc-denoisedSmoothed_bold.dtseries.nii'];
elseif string(task)=="rs2"
        fp=[parentfp '/' subj '_' sesh '_task-rs_acq-mb_dir-pe0_run-0_space-fsLR_den-91k_desc-denoisedSmoothed_bold.dtseries.nii'];
else
        fp=[parentfp '/' subj '_' sesh '_task-' task '_acq-mb_dir-pe0_run-0_space-fsLR_den-91k_desc-denoisedSmoothed_bold.dtseries.nii'];
end
% read it in
cifti_real = cifti_read(fp);
cifti2 = cifti_real;

% load simulated data
simFP=['/scratch/users/apines/ciftiout_Sym_' subj '_' sesh '_' task '.dtseries.nii'];
cifti_Sim = cifti_read(simFP);

% scale according to SD
% get maxes and mins to match real data
fmrStd=mean(std(cifti_real.cdata));
simStd=mean(std(cifti_Sim.cdata));

% calculate the factor to rescale the timeseries at
ReScaleFactor=fmrStd/simStd;
timecourse_sim=cifti_Sim.cdata.*ReScaleFactor;

% correct length
sizeTS=size(cifti_Sim.cdata);
cifti2.diminfo{2}.length=sizeTS(2);

cifti2.cdata = timecourse_sim;
cifti_write(cifti2,['/scratch/users/apines/ciftiout_Sym_' subj '_' sesh '_' task '.dtseries.nii']);
