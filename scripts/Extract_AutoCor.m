function Extract_AutoCor(subj,sesh,task)

% set parent directory
parentfp=['/scratch/groups/leanew1/xcpd_outP50_36p_bp/xcp_d/' subj '/' sesh '/func'];

% define some paths 
Paths{1} = '/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti/tool_folder';
addpath(genpath(Paths{1}))

% get autocorrelation FWHM for time series
C=ft_read_cifti_mod(['/scratch/groups/leanew1/xcpd_outP50_36p_bp/bv_concat/' Subject '-concat.dtseries.nii'])
C=ft_read_cifti_mod([parentfp '/' subj '_' sesh '_task-' task '_acq-mb_dir-pe0_run-0_space-fsLR_den-91k_desc-denoised_bold.dtseries.nii']; 
timeseries=C.data;

% get DMN-masked time series 
% load in DMN
dmn=ft_read_cifti_mod('/oak/stanford/groups/leanew1/users/apines/data/RobustInitialization/SoftParcel.dscalar.nii');
% mask as 0.3 or greater
dmn.data(dmn.data<0.3)=0;
dmn.data(dmn.data>=0.3)=1;
% apply to time series
%%% START HERE

% Define the number of time points and the number of time series
[num_verts, num_tps] = size(timeseries);
% Preallocate an array to store FWHM values for each time series
fwhm_values = zeros(num_verts, 1);

% Compute the autocorrelation FWHM for each time series
for i = 1:num_verts
    i
    % Extract the current time series
    current_time_series = timeseries(i, :);

    % Compute the autocorrelation of the current time series
    [autocorr, lags] = xcorr(current_time_series,'normalized',10);

    interpolated_lags = linspace(min(lags), max(lags), 100); % Adjust the number of points as needed
    interpolated_autocorr = interp1(lags, autocorr, interpolated_lags);

    % Determine the indices where autocorrelation is greater than or equal to half of the maximum value
    half_max_value = .5;
    above_half_max = interpolated_autocorr >= half_max_value;

    % Find the first and last indices where autocorrelation is above half-maximum
    first_index = find(above_half_max, 1, 'first');
    last_index = find(above_half_max, 1, 'last');

    % Calculate the FWHM for the current time series
    fwhm = last_index - first_index + 1; % FWHM in terms of number of data points

    % Store the FWHM value in the array, convert to seconds
    fwhm_values(i) = fwhm*.355;
end

% insert DMN-only vertices into broader grayordinate vector
% make full-size grayordinate vector
full_grayord_vector=zeros(91282,1);
% insert DMN-only vertices into full grayordinate vector
full_grayord_vector(dmn.data==1)=fwhm_values;
% save out
InfoMap_Ci = read_cifti('/oak/stanford/groups/leanew1/users/apines/NeuroSynthMaps/dmn_smooth.dscalar.nii')
InfoMap_Ci.cdata=fwhm_values;
InfoMap_Ci.diminfo{2}.length=1
InfoMap_Ci.diminfo{2}.type='series'
InfoMap_Ci.diminfo{2}.seriesUnit='SECOND'
InfoMap_Ci.diminfo{2}.seriesStep=.71
InfoMap_Ci.diminfo{2}.seriesStart=0.00

write_cifti(InfoMap_Ci,['~/testACF.dtseries.nii']);


% FINISH HERE
% sum autocorrelation across grayordinates

% save out autocorrelation
