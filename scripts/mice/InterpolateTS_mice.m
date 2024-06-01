function InterpolateTS_mice(subj,run)
% this function is to interpolate the functional time series to within-segment and between-TR timepoints. Should be v stable in band-passed fMR signal
addpath(genpath('/oak/stanford/groups/leanew1/users/apines/libs/'))

% copying OpFl code for reading in
basefp='/scratch/users/apines/p50_mice/proc/20200228/'
% load in specified scan
if run==1
        fn = [basefp 'thy1gc6s_0p3mgkg_' subj '_preLSD0p3mgkg_1/masked_dff_Gro_Masked_Sml_BP_Smoothed_Sml.h5']
        if exist(fn)
                data=h5read(fn, '/processed_data');
                formask=h5read(fn, '/mask');
        else
                disp('no run found')
        end
end
%% add if/else
% post 1
if run==2
        fn = [basefp 'thy1gc6s_0p3mgkg_' subj '_postLSD0p3mgkg_0/masked_dff_Gro_Masked_Sml_BP_Smoothed_Sml.h5']
        if exist(fn)
                data=h5read(fn, '/processed_data');
                formask=h5read(fn, '/mask');
        else
                disp('no run found')
        end
end
% post 2
if run==3
        fn = [basefp 'thy1gc6s_0p3mgkg_' subj '_postLSD0p3mgkg_5/masked_dff_Gro_Masked_Sml_BP_Smoothed_Sml.h5']
        if exist(fn)
                data=h5read(fn, '/processed_data');
                formask=h5read(fn, '/mask');
        else
                disp('no run found')
        end
end
% post 3
if run==4
        fn = [basefp 'thy1gc6s_0p3mgkg_' subj '_postLSD0p3mgkg_10/masked_dff_Gro_Masked_Sml_BP_Smoothed_Sml.h5']
        if exist(fn)
                data=h5read(fn, '/processed_data');
                formask=h5read(fn, '/mask');
        else
                disp('no run found')
        end
end
% post 4
if run==5
        fn = [basefp 'thy1gc6s_0p3mgkg_' subj '_postLSD0p3mgkg_15/masked_dff_Gro_Masked_Sml_BP_Smoothed_Sml.h5']
        if exist(fn)
                data=h5read(fn, '/processed_data');
                formask=h5read(fn, '/mask');
        else
                disp('no run found')
        end
end
% post 5
if run==6
        fn = [basefp 'thy1gc6s_0p3mgkg_' subj '_postLSD0p3mgkg_20/masked_dff_Gro_Masked_Sml_BP_Smoothed_Sml.h5']
        if exist(fn)
                data=h5read(fn, '/processed_data');
                formask=h5read(fn, '/mask');
        else
                disp('no run found')
        end
end

% make boolean mask if needed (use dmn mask later for refining)
mask = cellfun(@(x) strcmp(x, 'TRUE'), formask);

% get sizes
sizefMR=size(data);
lengthTS=sizefMR(3);

% store original timepoints as sequence
originalTimepoints = 1:lengthTS; 
% note this is the same in the mouse version because there's no need for scrubbing
timepointsToInterpolate = originalTimepoints;
% Calculate halfway timepoints
halfwayTimepoints = timepointsToInterpolate(1:end-1) + 0.5;
% initialize interpolated data
InterpData=zeros(sizefMR(1),sizefMR(2),length(halfwayTimepoints));
% get creative to deal with 3d structure here... loop over each column and each row, run interp timeseries? would avoid unraveling and raveling again
xdim=sizefMR(1);
ydim=sizefMR(2);
% loop over every pixel to interpolate signal over the time series
for x=1:xdim
	for y=1:ydim
		InterpData(x,y,:)=interp1(originalTimepoints, squeeze(data(x,y,:)), halfwayTimepoints);
	end
end

% save it out
OfInterpFp=[basefp '/' subj '_run-' num2str(run) '_interp.mat'];
save(OfInterpFp,'InterpData');

% load in DMN mask to create interpolated time series equivalent to output from extract_relative angles
networks=load('/oak/stanford/groups/leanew1/users/apines/p50_mice/Mouse_DMN_components.mat');
Dnet=networks.components_matrix';
% deploy equivalent masking procedure so we can have a matching vector field and interpolated time series for making histograms
nets=Dnet;
% create pixel-wise network mask
DMN_bool=Dnet;
DMN_bool(DMN_bool>.6)=1;
DMN_bool(DMN_bool<.6)=0;
DMN_bool=logical(DMN_bool);
% initialize matrix for each pixel to saveout to scratch
faceMatrix=zeros(sum(sum(DMN_bool)));
% network of interest
net=Dnet;
% calculate gradient of DMN
[nGx, nGy]=imgradientxy(net);
% use DMN threshold to binarize the DMN gradient
nGx=nGx(DMN_bool);
nGy=nGy(DMN_bool);

% pull out DMN pixels from interpolated time series using equivalent masking protocol (equivalent to that used in Extract_RelativeAngles_mice.m)
% note -1 is because optical flow occurs between frames, rendering the time series shorter by 1 (2 frames would have 1 OF measurement)
InterpData_DMN=zeros(sum(sum(DMN_bool)),lengthTS-1);
for F=1:length(nGx);
	% ineffcient, but matches other script so potentially more clear/consistent
	for tp=1:(lengthTS-1)
		% pull data for this one timepoint
		curData=InterpData(:,:,tp);
		currDataDMN=curData(DMN_bool);
		currPixel=currDataDMN(F);		
		InterpData_DMN(F,tp)=currPixel;
	end
end
% saveout
ofp=[basefp '/' subj '_run-' num2str(run) '_masked_interp_DMN.csv'];
% saveout
csvwrite(ofp,InterpData_DMN)

