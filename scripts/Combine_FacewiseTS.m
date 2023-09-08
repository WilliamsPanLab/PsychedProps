function Combine_FacewiseTS(subj,sesh)
% load in interpolated facewise time series, combine with angular distance time series
childfp=['/scratch/users/apines/data/mdma/' subj '/' sesh ];
ofpl=[childfp '/' subj '_' sesh '_task-rs_p2mm_masked_interp_L_faces.csv'];
ofpr=[childfp '/' subj '_' sesh '_task-rs_p2mm_masked_interp_R_faces.csv'];
% load in interpTS
iTS_l=readtable(ofpl);
iTS_r=readtable(ofpr);
% read in angular distances
aTS_l=readtable([childfp '/' subj '_' sesh '_Prop_TS_dmn_L.csv']);
aTS_r=readtable([childfp '/' subj '_' sesh '_Prop_TS_dmn_R.csv']);
% get dimensions to combine on
spatialDim=size(iTS_l);
spatialDim=spatialDim(1);
temporalDim=size(iTS_l);
temporalDim=temporalDim(2);
% combine
FullMatrix_L=zeros(spatialDim,temporalDim,2);
FullMatrix_L(:,:,1)=table2array(iTS_l);
FullMatrix_L(:,:,2)=table2array(aTS_l);
% right hemi
spatialDim=size(iTS_r);
spatialDim=spatialDim(1);
temporalDim=size(iTS_r);
temporalDim=temporalDim(2);
FullMatrix_R=zeros(spatialDim,temporalDim,2);
FullMatrix_R(:,:,1)=table2array(iTS_r);
FullMatrix_R(:,:,2)=table2array(aTS_r);
% DMN mask
% load in DMN
% subject to save masking as prop TS and interp_faces
% saveout
ofl=[childfp '/' subj '_' sesh '_task-rs_p2mm_masked_Bold_and_Angles_L.mat'];
ofr=[childfp '/' subj '_' sesh '_task-rs_p2mm_masked_Bold_and_Angles_R.mat'];
save(ofl, 'FullMatrix_L');
save(ofr, 'FullMatrix_R');
