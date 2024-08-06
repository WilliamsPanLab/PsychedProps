function Extract_AutoCor_psil(subj,sesh,task)

% set parent directory
parentfp=['/scratch/users/apines/PsiloData/' subj '/' subj '_' sesh '/func/' subj '_' sesh '_' task '.dtseries.nii'];

% define some paths 
Paths{1} = '/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti/tool_folder';
addpath(genpath(Paths{1}))

% read in citi
C=ft_read_cifti_mod(parentfp);

% extract time series
timeseries=C.data;

% read in gordon parcellation
GP=ft_read_cifti_mod('~/null_lL_WG33/Gordon333_FreesurferSubcortical.32k_fs_LR.dlabel.nii');

% gordon parcel labels
% https://balsa.wustl.edu/file/JX5V

% get DMN-masked time series 
% load in DMN
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
% initialize autocor vector
ACDMN=[];

% for each parcel, get complexity of timeseries
for p=DMNParcels
	% within each parcel
	% Define the number of time points and the number of time series
	[num_verts, num_tps] = size(timeseries(GP.data==p,:));
	% Preallocate an array to store FWHM values for each time series
	fwhm_values = zeros(num_verts, 1);
	% vertex indices
	inds=(GP.data==p);
	% Compute the autocorrelation FWHM for each time series
	for i = 1:num_verts
	    i;
	    % get this particular vertex
	    ind=inds(i);
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

	    % Store the FWHM value in t`he array, convert to seconds, divided by approx. 5 to account for interpolation above
	    fwhm_values(i) = (fwhm*.71)/5;
	end
	% get average autocor in this parcel
	avAC=mean(fwhm_values);
	ACDMN = [ACDMN avAC];
end

% save out autocor for dmn
avAC=mean(ACDMN);
T=table(avAC,'RowNames',"Row1");
outFP=['/scratch/users/apines/data/psil/' subj '/' sesh];
writetable(T,[outFP '/' subj '_' sesh '_' task '_AutoCor_gro.csv'],'WriteRowNames',true)


