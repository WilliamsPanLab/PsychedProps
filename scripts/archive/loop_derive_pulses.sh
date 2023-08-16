#ml python/3
for subj in sub-MDMA001 sub-MDMA002 sub-MDMA003 sub-MDMA005 sub-MDMA006 sub-MDMA007 sub-MDMA008 sub-MDMA009 sub-MDMA010 sub-MDMA011 sub-MDMA012 sub-MDMA013 sub-MDMA014 sub-MDMA015 sub-MDMA016 sub-MDMA017
do	
	for sesh in ses-00 ses-01 ses-02 ses-03
	do
		#/oak/stanford/groups/leanew1/users/apines/scripts/OpFl_CDys/scripts/sbatch_TheWorksReport.sh $subj $sesh
		#python3 derive_powers.py $subj $sesh
		sbatch /oak/stanford/groups/leanew1/users/apines/scripts/OpFl_CDys/scripts/fs_5/sbatch_OpFlow_MDMA.sh $subj $sesh
	done
done
