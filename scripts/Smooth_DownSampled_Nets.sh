subj=$1
# note: input "Group" if performing at group instead of subject-level
ml biology
ml workbench
# set directory where networks are at (and output to be)
netdir=/oak/stanford/groups/leanew1/users/apines/data/Atlas_Visualize

# for each network
for i in {1..4};
	do
	LinNet=${netdir}/${subj}_lh_Network_${i}_3k.func.gii
	LoutNet=${netdir}/${subj}_lh_Smooth_${i}_3k.func.gii
	RinNet=${netdir}/${subj}_rh_Network_${i}_3k.func.gii
	RoutNet=${netdir}/${subj}_rh_Smooth_${i}_3k.func.gii

	# use connectome workbench to smooth the networks 
	wb_command -metric-smoothing /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fsaverage4_std_sphere.L.3k_fsavg_L.surf.gii ${LinNet} 10 ${LoutNet}
	wb_command -metric-smoothing /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fsaverage4_std_sphere.R.3k_fsavg_R.surf.gii ${RinNet} 10 ${RoutNet}

done
