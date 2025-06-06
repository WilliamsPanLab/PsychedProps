function Extract_DMNSeg_lsd(subj,sesh,task)
% set parent directory
parentfp=['/scratch/users/apines/LSD_ICL/rest_proc/' subj '/' subj '_' sesh '_' task '_filt.dtseries.nii'];

% define some paths 
Paths{1} = '/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti/tool_folder';
addpath(genpath(Paths{1}))

% read in citi
C=ft_read_cifti_mod(parentfp);

% extract time series
C_timeseries=C.data;

% use already-calculated motion mask to temporally censor the full-res cifti
% load in temporal mask
childfp=['/scratch/users/apines/LSD_ICL/rest_proc/' subj];
tmaskfp=[childfp '/' subj '_' sesh '_task-' task '_AllSegments.txt'];
tmask=load(tmaskfp);
% Loop through each row in Absolut
for row = 1:size(tmask, 1)
        if tmask(row, 3) == 1
                % Extract the start and end values from the current row
                startValue = tmask(row, 1);
                endValue = tmask(row, 2);
                % change TRwise_mask_cont to 1 where this sequence of continuous good TRs occurs
                TRwise_mask_cont(startValue:endValue)=1;
        else
        end
end

% get correlation matrix of full time series... just cortex
C_timeseries=C_timeseries(1:59412,logical(TRwise_mask_cont));
% note cortical indices align even if subcort doesn't

% load in DMN
DMN=ft_read_cifti_mod('/oak/stanford/groups/leanew1/users/apines/maps/Network1_fslr.dscalar.nii');
DMNInds=find(DMN.data>.3);

% get correlation matrix of full time series... just cortex
C_timeseries=C_timeseries(1:59412,:);

% get non-dmn indices
nonDMNinds=setdiff(1:59412,DMNInds);
DMN_mat=C_timeseries(DMNInds,:);
nonDMN_mat=C_timeseries(nonDMNinds,:);
% avoid making full correlation matrix due to memory demands
correlation_matrix = 1 - pdist2(DMN_mat, nonDMN_mat, 'correlation');  
BWFC=mean(mean(correlation_matrix,'omitnan'),'omitnan');

% save out
T=table(BWFC,'RowNames',"Row1");
outFP=['/scratch/users/apines/LSD_ICL/rest_proc/' subj];
writetable(T,[outFP '/' subj '_' sesh '_' task '_DMNSeg.csv'],'WriteRowNames',true)


