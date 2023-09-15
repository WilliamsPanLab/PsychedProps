#!/bin/bash
#
#SBATCH --job-name=OpFl_Sim
#SBATCH --time=40:00:00
#SBATCH -n 9
#SBATCH --mem=45G
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

# downsample it 
/oak/stanford/groups/leanew1/users/apines/scripts/OpFl_CDys/scripts/DS_surf_ts_Simulated_fs5.sh $seed

# Calculate Optical Flow
matlab -nodisplay -r "OpFl_simulated_fs5('$seed')"

# simulate streamlines - left
matlab -nodisplay -r "OpFlStreamlines_Sim_L('$seed')"
matlab -nodisplay -r "OpFlStreamlines_Sim_R('$seed')"

# remove interim files
AgTS=/scratch/users/apines/ciftiout_Sym_${seed}.dtseries.nii
LeftHemi=/scratch/users/apines/Sim_L_AggTS_${seed}.func.gii
RightHemi=/scratch/users/apines/Sim_R_AggTS_${seed}.func.gii
LeftHemi_10=/scratch/users/apines/Sim_L_AggTS_${seed}_10k.func.gii
RightHemi_10=/scratch/users/apines/Sim_R_AggTS_${seed}_10k.func.gii
LeftHemi_10_mgh=/scratch/users/apines/Sim_L_AggTS_${seed}_10k.mgh
RightHemi_10_mgh=/scratch/users/apines/Sim_R_AggTS_${seed}_10k.mgh
# remove simulated TS
rm ${AgTS}
# remove sep. hemis
rm ${LeftHemi}
rm ${RightHemi}
# remove sep. hemis resampled
rm ${LeftHemi_10}
rm ${RightHemi_10}
# remove sep. hemis resampled .mgh
rm ${LeftHemi_10_mgh}
rm ${RightHemi_10_mgh}
