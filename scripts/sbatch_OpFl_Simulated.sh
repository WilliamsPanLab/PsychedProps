#!/bin/bash
#
#SBATCH --job-name=OpFl_Sim
#SBATCH --time=40:00:00
#SBATCH -n 4
#SBATCH --mem=30G
#SBATCH -p normal,leanew1  # Queue names you can submit to
# Outputs ----------------------------------
#SBATCH --mail-user=apines@stanford.edu
#SBATCH --mail-type=ALL
# ------------------------------------------


############    Buffaflow   ###################
##             _.-````'-,,
##   _,.,_ ,-'`           ``''-.,
## /)     (\                   '``-.
##((      ) )                ->    `\
## \)    (_/               ->       )\
##  |       /)           '    ,'    / \
##  `\    ^'            '     (    /  ))
##    |      _/\ ,     /    ,,`\   (  "`
##     \Y,   |  \  \  /____| / \_ \##
# This script is to prep the environment of a slurm node, 
##       `)_/    \  \  )    (->  (->
# Launch optical flow, caclulate angular distances, generate figures
##                \( \(     |/   |/
# and delete interim files upon completion
##    mic & dwb  /_(/_(    /_(  /_(

# random seed as argument
seed=$1

# will need matlab
module load matlab
# and freesurfer
module load biology
module load freesurfer/7.3.2
# and workbench
module load workbench
module load contribs poldrack anaconda/5.0.0-py36

# simulate data
matlab -nodisplay -r "SimulateFMR('$seed')"

# smooth it 
/oak/stanford/groups/leanew1/users/apines/scripts/OpFl_CDys/scripts/Smooth_Sim.sh $seed

# make sure we are all caught up
sleep 8

# convert to same value range
matlab -nodisplay -r "Scale_Simulated('$seed')"

# downsample it 
/oak/stanford/groups/leanew1/users/apines/scripts/OpFl_CDys/scripts/DS_surf_ts_Simulated.sh $seed

# Calculate Optical Flow
matlab -nodisplay -r "OpFl_simulated('$seed')"

# simulate streamlines - left
matlab -nodisplay -r "OpFlStreamlines_Sim_L('$seed')"
#matlab -nodisplay -r "OpFlStreamlines_Sim_R('$seed')"

# remove interim files
AgTS=/scratch/users/apines/ciftiout_Sym_${seed}.dtseries.nii
LeftHemi=/scratch/users/apines/Sim_L_AggTS_${seed}.func.gii
RightHemi=/scratch/users/apines/Sim_R_AggTS_${seed}.func.gii
LeftHemi_3=/scratch/users/apines/Sim_L_AggTS_${seed}_3k.func.gii
RightHemi_3=/scratch/users/apines/Sim_R_AggTS_${seed}_3k.func.gii
LeftHemi_3_mgh=/scratch/users/apines/Sim_L_AggTS_${seed}_3k.mgh
RightHemi_3_mgh=/scratch/users/apines/Sim_R_AggTS_${seed}_3k.mgh
# remove simulated TS
rm ${AgTS}
# remove sep. hemis
rm ${LeftHemi}
rm ${RightHemi}
# remove sep. hemis resampled
rm ${LeftHemi_3}
rm ${RightHemi_3}
# remove sep. hemis resampled .mgh
rm ${LeftHemi_3_mgh}
rm ${RightHemi_3_mgh}
