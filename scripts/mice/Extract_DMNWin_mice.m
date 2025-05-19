function Extract_DMNSeg_mice(subj,run)
% set parent directoryP load in data: pre LSD
basefp='/oak/stanford/groups/leanew1/users/apines/p50_mice/proc2/proc/20200228/'
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
% and an "over time" version of the mask to apply to the time series
%fullMask = repmat(Mask, [1, 1, lenOpFl]);

% pull out non-DMN pixels
NotDnet=networks.mask';
NotDnet_Bool=NotDnet;
% binarization required because of downsampling
NotDnet_Bool(NotDnet_Bool>.5)=1;
NotDnet_Bool(NotDnet_Bool<.5)=0;
% remove DMN to get non-dmn brain regions
NotDnet_Bool(DMN_bool==1)=0;
% add "over time" version
NotDMN_Mask=logical(NotDnet_Bool);
%full_NDMN_Mask = repmat(NotDMN_Mask, [1, 1, lenOpFl]);

% initialize DMN time series
DMN_mat = zeros(sum(Mask(:)), lenOpFl);
% use for loop to straightforwardly extract time series in pixels of interest
for t = 1:lenOpFl
    ts_slice = C_timeseries(:, :, t);
    DMN_mat(:, t) = ts_slice(Mask);
end

% initialize non DMN time series
nonDMN_mat = zeros(sum(NotDMN_Mask(:)), lenOpFl);
% use for loop to straightforwardly extract time series in pixels of interest
for t = 1:lenOpFl
    ts_slice = C_timeseries(:, :, t);
    nonDMN_mat(:, t) = ts_slice(NotDMN_Mask);
end

% avoid making full correlation matrix due to memory demands
correlation_matrix = 1 - pdist2(DMN_mat, DMN_mat, 'correlation');
% get average
WINFC=mean(mean(correlation_matrix,'omitnan'),'omitnan');

% save out dmn segregation
avDMNWIN=WINFC;

T=table(avDMNWIN,'RowNames',"Row1");
% makedir
outFP=['/scratch/users/apines/data/mouse/' subj '/'];
% makedir
writetable(T,[outFP '/' subj '_' num2str(run) '_avDMNWin.csv'],'WriteRowNames',true)


