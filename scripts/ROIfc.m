function ROIfc(subject,session)
	% derive full-brain fc from region of interest
	% subject: subject id
	% session: session id

% add nec. paths
addpath(genpath('/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti'));
% load in parcellated time series
childfp=['/scratch/users/apines/data/mdma/' subj '/' sesh '/'];
SubcortTS1=[childfp subj '_' sesh '_rs_SubCortROIS.txt'];
SubcortTS2=[childfp subj '_' sesh '_rs2_SubCortROIS.txt'];
SubcortTS1=readtable(SubcortTS1);
SubcortTS2=readtable(SubcortTS2);
% convert to arrays
SubcortTS1=table2array(SubcortTS1);
SubcortTS2=table2array(SubcortTS2);
% concatenate subcortical time series
SubcortTS=horzcat(SubcortTS1,SubcortTS2);
% load in time series
vw_ts_l_p=['/scratch/users/apines/data/mdma/' subj '/' sesh '/' subj '_' sesh '_L_AggTS_3k.mgh'];
vw_ts_r_p=['/scratch/users/apines/data/mdma/' subj '/' sesh '/' subj '_' sesh '_R_AggTS_3k.mgh'];
vw_ts_l=MRIread(vw_ts_l_p);
vw_ts_r=MRIread(vw_ts_r_p);
vw_ts_l=vw_ts_l.vol;
vw_ts_r=vw_ts_r.vol;
% load in medial wall mask: can include SNR mask here as well
surfML = '/oak/stanford/groups/leanew1/users/apines/surf/lh.Medial_wall.label';
mwIndVec_l = read_medial_wall_label(surfML);
surfMR = '/oak/stanford/groups/leanew1/users/apines/surf/rh.Medial_wall.label';
mwIndVec_r = read_medial_wall_label(surfMR);
% 6 additional vertices per hemisphere are decimated in the downsample
zerosL=find(all(vw_ts_l==0,4));
zerosR=find(all(vw_ts_r==0,4));
unionMask_L=union(zerosL,mwIndVec_l);
unionMask_R=union(zerosR,mwIndVec_r);
vw_ts_l_masked=vw_ts_l(1,setdiff([1:2562],unionMask_L),1,:);
vw_ts_r_masked=vw_ts_r(1,setdiff([1:2562],unionMask_R),1,:);
% squeeze to reduce extra dimensions
vw_ts_l_masked=squeeze(vw_ts_l_masked);
vw_ts_r_masked=squeeze(vw_ts_r_masked);

% get grayordinate-wise correlations with ROI timeseries
% left hemisphere
leftCorr=corr(SubcortTS(34,:)',vw_ts_l_masked')
% right hemisphere
rightCorr=corr(SubcortTS(9,:)',vw_ts_r_masked');
% save out FC maps as csv
csvwrite([childfp subj '_' sesh '_L_ROIfc.csv'],leftCorr);
csvwrite([childfp subj '_' sesh '_R_ROIfc.csv'],rightCorr);

end
