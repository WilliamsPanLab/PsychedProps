function us = OpFl_mdma_fs5_AmpSD_TS(subj,sesh)
% add paths
addpath(genpath('/oak/stanford/groups/leanew1/users/apines/libs/'))

% set filepath
childfp=['/scratch/users/apines/data/mdma/' subj '/' sesh ];
% load in processed opfl time series
OFtsL=[childfp '/OpFl_timeseries_L_fs5.mat'];
OFtsR=[childfp '/OpFl_timeseries_R_fs5.mat'];

% average over amplitude time series

% SD over amplitude: sliding window?

% generate grayplots of TD distance per-vert x PG positioning of vertex?

% add extra 0 as first timepoint: that way vectors represent movement from previous image

% saveout, not with this name and perhaps as pngs

save([childfp '/' subj '_' sesh '_OpFl_rs_fs5.mat'],'us')
save(outFP_L,'LvertTS_numfs5','-v7.3')
