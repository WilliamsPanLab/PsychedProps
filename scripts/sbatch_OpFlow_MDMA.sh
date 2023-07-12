#!/bin/bash
#
#SBATCH --job-name=OpFl
#SBATCH --time=5:00:00
#SBATCH -n 1
#SBATCH --mem=85G
#SBATCH -p leanew1 # Queue names you can submit to
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
rsOut=${childfp}/${subj}_${sesh}_NetDist_rs.mat

# Convert networks to dscalar
matlab -nodisplay -r "mat_to_dscalar('$subj')"

# Downsample networks
/oak/stanford/groups/leanew1/users/apines/scripts/FunctionalImages/Networks/DS_surf_Networks.sh $1

# Convert them back to mats for ang distance calc (in matlab)
matlab -nodisplay -r "Netgiis_2_mat('$subj')"

# make output directory outside scratch
mkdir /oak/stanford/groups/leanew1/users/apines/OpFlAngDs/mdma/${subj} 

matlab -nodisplay -r "Extract_RelativeAngles('$subj','$rsIn')"

