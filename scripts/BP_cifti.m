function BP_cifti(subj,sesh,task)
% addpaths
ToolFolder='/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti/tool_folder';
addpath(genpath(ToolFolder));

% Load CIFTI
parentfp=['/scratch/users/apines/LSD_ICL/rest_proc/' subj ]; 
infile = [parentfp '/' subj '_' sesh '_' task '.dtseries.nii'];
outname = [parentfp '/' subj '_' sesh '_' task '_filt.dtseries.nii'];
TR = 2;
Fs = 1 / TR;
N = 2;            % 2nd-order butterworth to match
% filter range matched
low_cutoff = 0.009;
high_cutoff = 0.08;
% load in
ts = read_cifti(infile);
data = ts.cdata;
n_time = size(data, 2);

% crete filter
[b, a] = butter(N, [low_cutoff high_cutoff] / (Fs / 2), 'bandpass');

% Filter grayordinates
data_filt = data;   % copy to preserve structure
for g = 1:size(data,1)
    signal = data(g, :);
    if any(isnan(signal)) || all(signal==0)
        continue  % skip bad rows
    end
    signal_filt = filtfilt(b, a, double(signal));
    data_filt(g, :) = signal_filt;
end

% Replace dtseries with filtered version
ts.cdata = data_filt;

% Save to new CIFTI file
write_cifti(ts,outname);

