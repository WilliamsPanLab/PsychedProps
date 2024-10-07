function vis_fs4PG(subj,sesh,fn)
% use same masking procedure as fc calc to visualize dmap output
% add filepaths
addpath(genpath('/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti'))
%fs4 fp
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
% load in PG: fill it in @ (setdiff[1:2562],unionMask_L) 
PG=load(['/scratch/users/apines/data/mdma/' subj '/' sesh '/' subj '_' sesh '_vertexwise_PG1.mat']);
PG=PG.pg;
PG_L=zeros(1,2562);
PG_L(setdiff([1:2562],unionMask_L))=PG(1:length(setdiff([1:2562],unionMask_L)));
PG_R=zeros(1,2562);
PG_R(setdiff([1:2562],unionMask_R))=PG((length(setdiff([1:2562],unionMask_L))+1):end);
Vis_Vertvec(PG_L,PG_R,fn)
