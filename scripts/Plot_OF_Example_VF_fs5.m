addpath(genpath('/oak/stanford/groups/leanew1/users/apines/libs/'))
% load in angles
parentFP=['/scratch/users/apines/data/mdma/' subj '/' sesh]
fpl=[parentFP '/' subj '_' sesh '_task-' task '_p2mm_masked_Vert_Angles_L.mat'];
fpr=[parentFP '/' subj '_' sesh '_task-' task '_p2mm_masked_Vert_Angles_R.mat'];
vecsL=load(fpl).vertWise_Vecs_l;
vecsR=load(fpr).vertWise_Vecs_r;
% run vis_vec_only_fs5 for a few timepoints
Vis_Vec_only_fs5(vecsL(:,1,:),vecsR(:,1,:),['~/' subj '_' sesh '_' task '_Mags1.png']);
Vis_Vec_only_fs5(vecsL(:,25,:),vecsR(:,25,:),['~/' subj '_' sesh '_' task '_Mags25.png']);
Vis_Vec_only_fs5(vecsL(:,75,:),vecsR(:,75,:),['~/' subj '_' sesh '_' task '_Mags75.png']);
Vis_Vec_only_fs5(vecsL(:,100,:),vecsR(:,100,:),['~/' subj '_' sesh '_' task '_Mags100.png']);
