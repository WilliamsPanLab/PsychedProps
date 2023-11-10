#!/bin/bash
#
#SBATCH --job-name=OpFl
#SBATCH --time=40:00:00
#SBATCH -n 4
#SBATCH --mem=20G
#SBATCH -p leanew1  # Queue names you can submit to
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

#############################
#### module III: Calc. Angles
#############################
echo "ΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔ"
echo "Starting module III: Angular distance calculation"
echo "ΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔ"
# make output directory outside scratch
mkdir /oak/stanford/groups/leanew1/users/apines/OpFlAngDs/mdma/${subj} 

# extract relative angles
# group
matlab -nodisplay -r "Extract_RelativeAngles('$subj','$sesh','rs1')"
matlab -nodisplay -r "Extract_RelativeAngles('$subj','$sesh','rs2')"
matlab -nodisplay -r "Extract_RelativeAngles('$subj','$sesh','emotion')"
matlab -nodisplay -r "Extract_RelativeAngles('$subj','$sesh','gambling')"
matlab -nodisplay -r "Extract_RelativeAngles('$subj','$sesh','wm')"

#################
echo "OpFl complete"
