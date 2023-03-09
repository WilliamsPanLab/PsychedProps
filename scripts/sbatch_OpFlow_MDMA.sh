#!/bin/bash
#
#SBATCH --job-name=OpFl
#SBATCH --time=5:00:00
#SBATCH -n 1
#SBATCH --mem=85G
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
# and workbench
module load workbench

# subject name is input argument
subj=$1
# sesh is input 2
sesh=$2

# Downsample the data 
/oak/stanford/groups/leanew1/users/apines/scripts/OpFl_CDys/scripts/DS_surf_ts_mdma.sh $1 $2

# cd to scripts directory
cd /oak/stanford/groups/leanew1/users/apines/scripts/OpFl_CDys/scripts

# Calculate Optical Flow
matlab -nodisplay -r "OpFl_mdma('$subj','$sesh')"

# RS filepaths
childfp=/scratch/users/apines/data/mdma/${subj}/${sesh}
rsIn=${childfp}/${subj}_${sesh}_OpFl_rs.mat
rsOut=${childfp}/${subj}_${sesh}_PGGDist_rs.mat

# Calculate Angular Distance: rs
/oak/stanford/groups/leanew1/users/apines/scripts/OpFl_CDys/scripts/run_PGG_AngDistCalc4_CompVer_PL.sh /share/software/user/restricted/matlab/R2018a/ $rsIn $rsOut

# make output directory outside scratch
 mkdir /oak/stanford/groups/leanew1/users/apines/OpFlAngDs/mdma/${subj} 

#./run_Extract_BUTD_ResultantVecs_PL.sh /share/software/user/restricted/matlab/R2018a/ $rsIn $rsOut /oak/stanford/groups/leanew1/users/apines/OpFlAngDs/mdma/${subj}/${subj}_${sesh}_rs_BUTD_L.mat /oak/stanford/groups/leanew1/users/apines/OpFlAngDs/mdma/${subj}/${subj}_${sesh}_rs_BUTD_R.mat

# convert to R format
matlab -nodisplay -r "BUTD_to_Rformat_MDMA('$subj','$sesh')"

./run_Extract_BUTD_ResultantVecs_Gran.sh /share/software/user/restricted/matlab/R2018a/ $rsIn $childfp/OpFl_timeseries_L.mat $childfp/OpFl_timeseries_R.mat
