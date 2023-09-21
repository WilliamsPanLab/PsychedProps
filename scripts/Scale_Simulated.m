function Scale_Simulated(seed)
% scale simulated data to be same magnitude as real data (range of signal)
addpath(genpath('/oak/stanford/groups/leanew1/users/apines/libs/'))
% load real data
cifti_real = cifti_read('~/sub-MDMA006_ses-00_task-rs_acq-mb_dir-pe0_run-0_space-fsLR_den-91k_desc-denoisedSmoothed_bold.dtseries.nii');
% * 2 becaue each subj has two rs scans per sesh
timelength = cifti_real.diminfo{2}.length*2;
cifti2 = cifti_real;
cifti2.diminfo{2}.length=timelength;

% load simulated data
simFP=['/scratch/users/apines/ciftiout_Sym_' seed '.dtseries.nii'];
cifti_Sim = cifti_read(simFP);

% scale according to SD
% get maxes and mins to match real data
fmrStd=mean(std(cifti_real.cdata));
simStd=mean(std(cifti_Sim.cdata));

% calculate the factor to rescale the timeseries at
ReScaleFactor=fmrStd/simStd;
timecourse_sim=cifti_Sim.cdata.*ReScaleFactor;

cifti2.cdata = timecourse_sim;
cifti_write(cifti2,['/scratch/users/apines/ciftiout_Sym_' seed '.dtseries.nii']);
