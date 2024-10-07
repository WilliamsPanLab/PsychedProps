# and freesurfer
module load biology
module load freesurfer/7.3.2
# and workbench
module load workbench
# subject name is input argument
subj=$1
# sesh is input 2
sesh=$2

xcpd_outdir=/scratch/groups/leanew1/xcpd_outP50_36p_bp/xcp_d/${subj}/${sesh}/func/
rs1xcpSubcort_fp=${xcpd_outdir}${subj}_${sesh}_task-rs_acq-mb_dir-pe0_run-0_space-fsLR_den-91k_desc-smooth_alff.dscalar.nii
rs2xcpSubcort_fp=${xcpd_outdir}${subj}_${sesh}_task-rs_acq-mb_dir-pe0_run-1_space-fsLR_den-91k_desc-smooth_alff.dscalar.nii
childfp=/oak/stanford/groups/leanew1/users/apines/data/p50/${subj}/${sesh}/
alff_rs1xcpSubcort_fp=${xcpd_outdir}${subj}_${sesh}_task-rs_acq-mb_dir-pe0_run-0_space-fsLR_den-91k_desc-smooth_alff.dscalar.nii
alff_rs2xcpSubcort_fp=${xcpd_outdir}${subj}_${sesh}_task-rs_acq-mb_dir-pe1_run-0_space-fsLR_den-91k_desc-smooth_alff.dscalar.nii 
alff_p_SubcortTS1=${childfp}${subj}_${sesh}_rs_SubCortROIS_alff.ptseries.nii
alff_p_SubcortTS2=${childfp}${subj}_${sesh}_rs2_SubCortROIS_alff.ptseries.nii
alff_p_Subcort=${childfp}${subj}_${sesh}_rs_SubCortROIS_alff.ptseries.nii
alff_p_Subcort=${childfp}${subj}_${sesh}_rs2_SubCortROIS_alff.ptseries.nii
alff_SubcortTS1=${childfp}${subj}_${sesh}_rs_SubCortROIS_alff.txt
alff_SubcortTS2=${childfp}${subj}_${sesh}_rs2_SubCortROIS_alff.txt
scale3fp=/oak/stanford/groups/leanew1/users/apines/maps/Tian_Subcortex_S3_3T_32k.dlabel.nii
# convert alff to ptseries
wb_command -cifti-parcellate ${alff_rs1xcpSubcort_fp} ${scale3fp} COLUMN ${alff_p_SubcortTS1}
wb_command -cifti-parcellate ${alff_rs2xcpSubcort_fp} ${scale3fp} COLUMN ${alff_p_SubcortTS2}
# convert to text
wb_command -cifti-convert -to-text $alff_p_SubcortTS1 $alff_SubcortTS1
wb_command -cifti-convert -to-text $alff_p_SubcortTS2 $alff_SubcortTS2
