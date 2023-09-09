addpath(genpath('/oak/stanford/groups/leanew1/users/apines/libs/'))
%% adding realistic power spectra, adapted code from Tim Laumann's method to simulated power spectra from real data
% see Laumann et al. 2017 cer cortex paper
TRin = 1;  % randomly picked
TRout = 1;  % randomly picked

% run pwelch loop over regions
cifti_real = cifti_read('~/sub-MDMA006_ses-00_task-rs_acq-mb_dir-pe0_run-0_space-fsLR_den-91k_desc-denoisedSmoothed_bold.dtseries.nii');
% * 2 becaue each subj has two rs scans per sesh
timelength = cifti_real.diminfo{2}.length*2;
cifti2 = cifti_real;
cifti2.diminfo{2}.length=timelength;
cifti2.cdata =randn(timelength,91282);

parcelPS = [];
for reg = 1:1:91282
    reg
    [parcelPS(reg,:),freq] = pwelch(cifti_real.cdata(reg,:),[],[],[],1/TRin);
end

ps_mean = mean(parcelPS);
P_target = ps_mean;

wn = cifti2.cdata;

%%Normalize to unit variance
wn_std = std(wn);
wn_mean = mean(wn);
wn_unit = (wn-repmat(wn_mean,timelength,1))./repmat(wn_std,timelength,1);

%%Power spectral density resampling
if mod(timelength,2)
    midpoint = ceil(timelength./2);
else
    midpoint = timelength./2+1;
end

% Find frequency indices to interpolate
if mod(length(P_target),2)
    P_target_length = ceil(length(P_target)/2)+1;
else
    P_target_length = length(P_target)/2+1;
end

count = 1;
for k = 0:timelength/2;
    Vq(count) = (k*TRin*length(P_target))/(TRout*(timelength))+1;
    count = count + 1;
end

% Interpolate 
resamp_P_target = interp1(P_target(1:P_target_length),Vq);

% Flip over nyquist
resamp_P_target = [resamp_P_target fliplr(resamp_P_target(2:(midpoint-1)))]';

if mod(timelength,2) == 1
    resamp_P_target = [0; resamp_P_target];
end

%%Enforce power spectra
fft_timecourse_sim = fft(wn_unit);
resamp_P_target_rep = repmat(resamp_P_target,[1 size(fft_timecourse_sim,2)]);
F_sim = fft_timecourse_sim.*sqrt(resamp_P_target_rep);

timecourse_sim = ifft(F_sim,'symmetric');
timecourse_sim = (timecourse_sim-repmat(mean(timecourse_sim),size(timecourse_sim,1),1))./repmat(std(timecourse_sim),size(timecourse_sim,1),1);
cifti2.cdata = timecourse_sim';
cifti_write(cifti2,'~/ciftiout_Sym.dtseries.nii');
