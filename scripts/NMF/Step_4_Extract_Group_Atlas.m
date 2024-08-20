% add paths
ToolFolder='/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti/tool_folder';
addpath(genpath(ToolFolder));

% set directories
Folder = '/oak/stanford/groups/leanew1/users/apines/data/'

%% Group atlas was the clustering results of 50 atlases during the initialization
GroupAtlasLoading_Mat = load([Folder '/RobustInitialization/init.mat']);
SubjectsFolder = '/oak/stanford/groups/leanew1/users/apines/fs5surf/'

% for surface data
surfML = [SubjectsFolder '/lh.Mask_SNR.label'];
surfMR = [SubjectsFolder '/rh.Mask_SNR.label'];

% reconcile with medial wall
mwIndVec_l = read_medial_wall_label(surfML);
Index_l = setdiff([1:10242], mwIndVec_l);
mwIndVec_r = read_medial_wall_label(surfMR);
Index_r = setdiff([1:10242], mwIndVec_r);

% get network loadings
initV = GroupAtlasLoading_Mat.initV;
initV_Max = max(initV);
% omit negligible loadings
trimInd = initV < .5;
initV(trimInd) = 0;
sbj_AtlasLoading_NoMedialWall = initV;
[~, sbj_AtlasLabel_NoMedialWall] = max(sbj_AtlasLoading_NoMedialWall, [], 2);

% initialize output labels/loadings
sbj_AtlasLabel_lh = zeros(1, 10242);
sbj_AtlasLoading_lh = zeros(4, 10242);
% populate
sbj_AtlasLabel_lh(Index_l) = sbj_AtlasLabel_NoMedialWall(1:length(Index_l));
sbj_AtlasLoading_lh(:, Index_l) = sbj_AtlasLoading_NoMedialWall(1:length(Index_l), :)';
sbj_AtlasLabel_rh = zeros(1, 10242);
sbj_AtlasLoading_rh = zeros(4, 10242);
sbj_AtlasLabel_rh(Index_r) = sbj_AtlasLabel_NoMedialWall(length(Index_l) + 1:end);
sbj_AtlasLoading_rh(:, Index_r) = sbj_AtlasLoading_NoMedialWall(length(Index_l) + 1:end, :)';
% save out
save([Folder '/SingleAtlas_Analysis/Group_AtlasLabel.mat'], 'sbj_AtlasLabel_lh', 'sbj_AtlasLabel_rh', 'sbj_AtlasLabel_NoMedialWall');
save([Folder '/SingleAtlas_Analysis/Group_AtlasLoading.mat'], 'sbj_AtlasLoading_lh', 'sbj_AtlasLoading_rh', 'sbj_AtlasLoading_NoMedialWall');

% make a visualization folder
VisualizeFolder = [Folder '/Atlas_Visualize'];
mkdir(VisualizeFolder);

% convert to func.gii's
GroupAtlasLoading_Mat = load([Folder '/SingleAtlas_Analysis/Group_AtlasLoading.mat']);
for i = 1:4
  % left hemi
  V_lh = gifti;
  V_lh.cdata = GroupAtlasLoading_Mat.sbj_AtlasLoading_lh(i, :)';
  V_lh_File = [VisualizeFolder '/Group_lh_Network_' num2str(i) '.func.gii'];
  save(V_lh, V_lh_File);
  % right hemi
  V_rh = gifti;
  V_rh.cdata = GroupAtlasLoading_Mat.sbj_AtlasLoading_rh(i, :)';
  V_rh_File = [VisualizeFolder '/Group_rh_Network_' num2str(i) '.func.gii'];
  save(V_rh, V_rh_File);
  % convert into cifti file
  cmd = ['wb_command -cifti-create-dense-scalar ' VisualizeFolder '/Group_AtlasLoading_Network_' num2str(i) ...
         '.dscalar.nii -left-metric ' V_lh_File ' -right-metric ' V_rh_File];
  system(cmd);
  pause(1);
end

