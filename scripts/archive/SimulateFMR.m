function SimulateFMR(subj,sesh,task)
addpath(genpath('/oak/stanford/groups/leanew1/users/apines/libs/'))
%% adding realistic power spectra, adapted code from Tim Laumann's method to simulated power spectra from real data
rng(42069)

% see Laumann et al. 2017 cer cortex paper
TRin = .71; 
TRout = .71;

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
% run pwelch loop over regions
timelength = cifti_real.diminfo{2}.length;
cifti2 = cifti_real;
cifti2.diminfo{2}.length=timelength;
cifti2.cdata =randn(timelength,91282);

parcelPS = [];
for reg = 1:1:91282
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

% mask by subject-specific TRs before printing out
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% load in continuous segment indices
childfp=['/scratch/users/apines/data/mdma/' subj '/' sesh];
CSIfp=[childfp '/' subj '_' sesh '_task-' task '_AllSegments.txt'];
CSI = importdata(CSIfp);
% extract boolean of TRs we want
booleanVector = zeros(1, max(CSI(:, 2)));

% Iterate through the rows of CSI
for i = 1:size(CSI, 1)
    if CSI(i, 3) == 1
        % Set the elements in the range between column 1 and column 2 to 1
        booleanVector(CSI(i, 1):CSI(i, 2)) = 1;
    end
end

% insert back into OG structure to saveout in equivalent structure
cifti2.cdata = timecourse_sim(logical(booleanVector),:)';
% update length of sim cifti
cifti2.diminfo{2}.length=(sum(logical(booleanVector)));
% saveout
cifti_write(cifti2,['/scratch/users/apines/ciftiout_Sym_' subj '_' sesh '_' task '.dtseries.nii']);



