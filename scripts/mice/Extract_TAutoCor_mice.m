function Extract_TAutoCor_mice(subj,run)
% set parent directoryP load in data: pre LSD
basefp='/scratch/users/apines/p50_mice/proc/20200228/'
% load in specified scan
if run==1
        fn = [basefp 'thy1gc6s_0p3mgkg_' subj '_preLSD0p3mgkg_1/masked_dff_Gro_Masked_Sml_BP_Smoothed_Sml.h5']
        if exist(fn)
                data=h5read(fn, '/processed_data');
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
        else
                disp('no run found')
        end
end
% post 2
if run==3
        fn = [basefp 'thy1gc6s_0p3mgkg_' subj '_postLSD0p3mgkg_5/masked_dff_Gro_Masked_Sml_BP_Smoothed_Sml.h5']
        if exist(fn)
                data=h5read(fn, '/processed_data');
        else
                disp('no run found')
        end
end
% post 3
if run==4
        fn = [basefp 'thy1gc6s_0p3mgkg_' subj '_postLSD0p3mgkg_10/masked_dff_Gro_Masked_Sml_BP_Smoothed_Sml.h5']
        if exist(fn)
                data=h5read(fn, '/processed_data');
        else
                disp('no run found')
        end
end
% post 4
if run==5
        fn = [basefp 'thy1gc6s_0p3mgkg_' subj '_postLSD0p3mgkg_15/masked_dff_Gro_Masked_Sml_BP_Smoothed_Sml.h5']
        if exist(fn)
                data=h5read(fn, '/processed_data');
        else
                disp('no run found')
        end
end
% post 5
if run==6
        fn = [basefp 'thy1gc6s_0p3mgkg_' subj '_postLSD0p3mgkg_20/masked_dff_Gro_Masked_Sml_BP_Smoothed_Sml.h5']
        if exist(fn)
                data=h5read(fn, '/processed_data');
        else
                disp('no run found')
        end
end

% pull in data as "C", refers to cifti in human scripts
C_timeseries=data;
lenOpFl=size(C_timeseries,3);

% load in mouse DMN
networks=load('/oak/stanford/groups/leanew1/users/apines/p50_mice/Mouse_DMN_components.mat');
Dnet=networks.components_matrix';
DMN_bool=Dnet;
DMN_bool(DMN_bool>.6)=1;
DMN_bool(DMN_bool<.6)=0;
Mask=logical(DMN_bool);

% initialize DMN time series
DMN_mat = zeros(sum(Mask(:)), lenOpFl);
% use for loop to straightforwardly extract time series in pixels of interest
for t = 1:lenOpFl
    ts_slice = C_timeseries(:, :, t);
    DMN_mat(:, t) = ts_slice(Mask);
end

% initialize an autocor array
ACarray=zeros(1,sum(Mask(:)));
% number of pixels
NumPix=sum(Mask(:));

% loop over each pixel within DMN
for p = 1:NumPix;
	% extract data
	dataInPix=DMN_mat(p,:);
	% OG time series (-1)
	OGTS=dataInPix(1:end-1);
        ShiftTS=dataInPix(2:end);
	% calculate autocor (t-1 vs. t)
	correlationOfInterest=corrcoef(OGTS,ShiftTS);
        % note this returns a correlation matrix: only interested in the between-timeseries measurement (off-diag)
        ACarray(p)=correlationOfInterest(1,2);
% end pixel loop
end


% get average
av_AutoCor=mean(ACarray);

% saveout
T=table(av_AutoCor,'RowNames',"Row1");
% makedir
outFP=['/scratch/users/apines/data/mouse/' subj '/'];
% makedir
writetable(T,[outFP '/' subj '_' num2str(run) '_av_AutoCor.csv'],'WriteRowNames',true)


