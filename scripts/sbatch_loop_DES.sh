# note: be in folder w/ subj names as subfolders to run (so * picks up correct names)
for d in `cat ~/connectomeSubjs.txt`; 
do sbatch /oak/stanford/groups/leanew1/users/apines/scripts/OpFl_CDys/scripts/sbatch_DS_DES_fs5.sh $d;
done
