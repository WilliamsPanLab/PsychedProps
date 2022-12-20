function fs4FC(subj,sesh)
% add filepaths
addpath(genpath('/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti'))
%fs4 fp
vw_ts_l_p=['/scratch/users/apines/data/mdma/' subj '/' sesh '/' subj '_' sesh '_L_AggTS_3k.mgh'];
vw_ts_r_p=['/scratch/users/apines/data/mdma/' subj '/' sesh '/' subj '_' sesh '_R_AggTS_3k.mgh'];
vw_ts_l=MRIread(vw_ts_l_p);
vw_ts_r=MRIread(vw_ts_r_p);
vw_ts_l=vw_ts_l.vol;
vw_ts_r=vw_ts_r.vol;
% apply SNR masks when developed
vw_ts_l_masked=vw_ts_l;
vw_ts_r_masked=vw_ts_r;
%vw_ts_l_masked=vw_ts_l(1,(logical(surfMask.l)),1,:);
%vw_ts_r_masked=vw_ts_r(1,(logical(surfMask.r)),1,:);
% stacking matrices so vertex number is doubled
vw_ts_both=[vw_ts_l_masked vw_ts_r_masked];
% remove extra dimensions
size(vw_ts_both)
vw_ts_both=squeeze(vw_ts_both)';
size(vw_ts_both)
% bigass connectivity matrix, takes 5 seconds or so to calc
ba_conmat=corrcoef(vw_ts_both);
size(ba_conmat);
% save these out to scratch
fcPath=['/scratch/users/apines/data/mdma/' subj '/' sesh '/' subj '_' sesh '_5k_FC.csv'];
csvwrite(fcPath,ba_conmat);
