ml biology
ml workbench
# set directory where networks are at (and output to be)
netdir=/oak/stanford/groups/leanew1/users/apines/data/RobustInitialization/
# use connectome workbench to smooth the networks 
wb_command -metric-smoothing /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fsaverage4_std_sphere.L.3k_fsavg_L.surf.gii ${netdir}group_L_AggNets_3k.func.gii 20 ${netdir}/group_L_AggNets_3kSmooth.func.gii
wb_command -metric-smoothing /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fsaverage4_std_sphere.R.3k_fsavg_R.surf.gii ${netdir}group_R_AggNets_3k.func.gii 20 ${netdir}/group_R_AggNets_3kSmooth.func.gii

