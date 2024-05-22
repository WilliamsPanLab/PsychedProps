function Extract_NGSC_mice(subj,sesh,run)
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
% and an "over time" version of the mask to apply to the time series
fullMask = repmat(Mask, [1, 1, lenOpFl]);

% get complexity of full time series
masked_time_series=C_timeseries(fullMask);
dmn_ts = reshape(masked_time_series, sum(Mask(:)), lenOpFl);

% visualization check
% initialize matrix to hold data
%projected_time_series = zeros(67, 70, lenOpFl);

% Fill in the projected time series using the mask
%for t = 1:lenOpFl
%    temp_frame = zeros(67, 70);
%    temp_frame(Mask) = dmn_ts(:, t);
%    projected_time_series(:, :, t) = temp_frame;
%end

% Now you can visualize the projected time series
% For example, visualize the first time point
%fig=figure;
%imagesc(projected_time_series(:, :, 1));
%hold on;
%colorbar;
%hold on;
%title('Projected Time Series at Time Point 1');
%caxis([min(min((projected_time_series(:,:,1)))), max(max(projected_time_series(:,:,1)))]);
%print(fig, '~/MaskedTimeseries', '-dpng', '-r600')

% explicitly using Josh's code, note scrubbing mask is TRwise_mask_cont
[~,~,~,~,EXPLAINED]=pca(dmn_ts);
EXPLAINED=EXPLAINED/100;
nGSC=-sum(EXPLAINED .* log(EXPLAINED))/log(length(EXPLAINED));
cxDMN = nGSC;

% save out normalized entropy for dmn
avComplexity=cxDMN;

T=table(avComplexity,'RowNames',"Row1");
% makedir
outFP=['/scratch/users/apines/data/mouse/' subj '/'];
system(['mkdir ' outFP]);
% makedir
outFP=['/scratch/users/apines/data/mouse/' subj '/' sesh];
system(['mkdir ' outFP]);
writetable(T,[outFP '/' subj '_' sesh '_' num2str(run) '_Complexity_gro.csv'],'WriteRowNames',true)


