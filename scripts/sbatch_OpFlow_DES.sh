#!/bin/bash
#
#SBATCH --job-name=OpFl
#SBATCH --time=2:30:00
#SBATCH -n 1
#SBATCH --mem=15G
#SBATCH -p leanew1,normal  # Queue names you can submit to
# Outputs ----------------------------------
#SBATCH --mail-user=no@stanford.edu
#SBATCH --mail-type=ALL
# ------------------------------------------

# will need matlab
module load matlab
# and freesurfer
module load biology
module load freesurfer/7.3.2

# subject name is input argument
subj=$1

# Downsample the data 
/oak/stanford/groups/leanew1/users/apines/scripts/OpFl_CDys/scripts/DS_surf_ts_DES.sh $1

# cd to scripts directory
cd /oak/stanford/groups/leanew1/users/apines/scripts/OpFl_CDys/scripts

# Center and bandpass the data
matlab -nodisplay -r "C_and_BP_DES('$subj')"

# Calculate Optical Flow
matlab -nodisplay -r "OpFl_DES('$subj')"

# WM filepaths
wmIn=/scratch/users/apines/data/DES/${subj}/OpFl_wm.mat
wmOut=/scratch/users/apines/data/DES/${subj}/PGGDist_wm.mat

# RS1 filepaths
rsIn=/scratch/users/apines/data/DES/${subj}/OpFl_rs.mat
rsOut=/scratch/users/apines/data/DES/${subj}/PGGDist_rs.mat

# RS2 filepaths
rs2In=/scratch/users/apines/data/DES/${subj}/OpFl_rs2.mat
rs2Out=/scratch/users/apines/data/DES/${subj}/PGGDist_rs2.mat

# Calculate Angular Distance: wm
/oak/stanford/groups/leanew1/users/apines/scripts/OpFl_CDys/scripts/run_PGG_AngDistCalc4_CompVer_PL.sh /share/software/user/restricted/matlab/R2018a/ $wmIn $wmOut

# Calculate Angular Distance: rs
/oak/stanford/groups/leanew1/users/apines/scripts/OpFl_CDys/scripts/run_PGG_AngDistCalc4_CompVer_PL.sh /share/software/user/restricted/matlab/R2018a/ $rsIn $rsOut

# rs2
/oak/stanford/groups/leanew1/users/apines/scripts/OpFl_CDys/scripts/run_PGG_AngDistCalc4_CompVer_PL.sh /share/software/user/restricted/matlab/R2018a/ $rs2In $rs2Out

# make output directory outside scratch
mkdir /oak/stanford/groups/leanew1/users/apines/OpFlAngDs/DES/${subj}

# Extract Directional information
./run_Extract_BUTD_ResultantVecs_PL.sh /share/software/user/restricted/matlab/R2018a/ $wmIn $wmOut /oak/stanford/groups/leanew1/users/apines/OpFlAngDs/DES/${subj}/wm_BUTD_L.mat /oak/stanford/groups/leanew1/users/apines/OpFlAngDs/DES/${subj}/wm_BUTD_R.mat

./run_Extract_BUTD_ResultantVecs_PL.sh /share/software/user/restricted/matlab/R2018a/ $rsIn $rsOut /oak/stanford/groups/leanew1/users/apines/OpFlAngDs/DES/${subj}/rs_BUTD_L.mat /oak/stanford/groups/leanew1/users/apines/OpFlAngDs/DES/${subj}/rs_BUTD_R.mat

./run_Extract_BUTD_ResultantVecs_PL.sh /share/software/user/restricted/matlab/R2018a/ $rs2In $rs2Out /oak/stanford/groups/leanew1/users/apines/OpFlAngDs/DES/${subj}/rs2_BUTD_L.mat /oak/stanford/groups/leanew1/users/apines/OpFlAngDs/DES/${subj}/rs2_BUTD_R.mat

# convert to R format
matlab -nodisplay -r "BUTD_to_Rformat_DES('$subj')"


# LEGACY COMMENTED OUT UNTIL DATA PROC'ED
# GNG filepaths
#gngIn=/scratch/users/apines/data/DES/${subj}/OpFl_GNG.mat
#gngOut=/scratch/users/apines/data/DES/${subj}/PGGDist_GNG.mat

# NCF filepaths
#ncfIn=/scratch/users/apines/data/DES/${subj}/OpFl_NCF.mat
#ncfOut=/scratch/users/apines/data/DES/${subj}/PGGDist_NCF.mat

# CON filepaths
#conIn=/scratch/users/apines/data/ispot/${subj}/OpFl_CON.mat
#conOut=/scratch/users/apines/data/ispot/${subj}/PGGDist_CON.mat

# Calculate Angular Distance: GNG
#/oak/stanford/groups/leanew1/users/apines/scripts/OpFl_CDys/scripts/run_PGG_AngDistCalc4_CompVer_PL.sh /share/software/user/restricted/matlab/R2018a/ $gngIn $gngOut

# Calculate Angular Distance: NCF
#/oak/stanford/groups/leanew1/users/apines/scripts/OpFl_CDys/scripts/run_PGG_AngDistCalc4_CompVer_PL.sh /share/software/user/restricted/matlab/R2018a/ $ncfIn $ncfOut

# Con faces
#/oak/stanford/groups/leanew1/users/apines/scripts/OpFl_CDys/scripts/run_PGG_AngDistCalc4_CompVer_PL.sh /share/software/user/restricted/matlab/R2018a/ $conIn $conOut

# Extract Directional information
#./run_Extract_BUTD_ResultantVecs_PL.sh /share/software/user/restricted/matlab/R2018a/ $gngIn $gngOut /oak/stanford/groups/leanew1/users/apines/OpFlAngDs/ispot/${subj}/gng_BUTD_L.mat /oak/stanford/groups/leanew1/users/apines/OpFlAngDs/ispot/${subj}/gng_BUTD_R.mat

#./run_Extract_BUTD_ResultantVecs_PL.sh /share/software/user/restricted/matlab/R2018a/ $ncfIn $ncfOut /oak/stanford/groups/leanew1/users/apines/OpFlAngDs/ispot/${subj}/ncf_BUTD_L.mat /oak/stanford/groups/leanew1/users/apines/OpFlAngDs/ispot/${subj}/ncf_BUTD_R.mat

#./run_Extract_BUTD_ResultantVecs_PL.sh /share/software/user/restricted/matlab/R2018a/ $conIn $conOut /oak/stanford/groups/leanew1/users/apines/OpFlAngDs/ispot/${subj}/con_BUTD_L.mat /oak/stanford/groups/leanew1/users/apines/OpFlAngDs/ispot/${subj}/con_BUTD_R.mat

# convert to R format
#matlab -nodisplay -r "BUTD_to_Rformat('$subj')"


