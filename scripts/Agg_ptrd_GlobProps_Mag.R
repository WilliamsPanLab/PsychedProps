### MANUALLY RECREATE SUBJ LIST
subjPrefix=rep('sub-PTRD0',7)
# no subj 4 or 15
subjSuffix=c('04','06','07','09','10','11','12')
subjList=paste0(subjPrefix,subjSuffix)

# subjvec : no 3
subjvec=c(4,6,7,9,10,11,12)

# create output vectors for group-level inference on delay in each ROI
pl_lAMY_r=NULL
pl_mAMY_r=NULL
pl_lAMY_l=NULL
pl_mAMY_l=NULL
ps_lAMY_r=NULL
ps_mAMY_r=NULL
ps_lAMY_l=NULL
ps_mAMY_l=NULL

# for each subj
for (s in 1:length(subjvec)){
	subj=subjList[s]
	# get output filepath
	childfp=paste0('/oak/stanford/groups/leanew1/users/apines/data/p50/',subj,'/') 
	Del_baselineFPs=list.files(path=(paste0(childfp,'ses-00','/')),pattern="MagMat.csv")
	Del_ps1FPs=list.files(path=(paste0(childfp,'ses-01','/')),pattern="MagMat.csv")
	# load in each non-drug csv
	for (b in 1:length(Del_baselineFPs)){
		filename=paste0(paste0(childfp,'ses-00','/'),Del_baselineFPs[b])
		file=read.table(filename,header=FALSE,sep=',')
		pl_lAMY_r=c(pl_lAMY_r,file[1,])
		pl_mAMY_r=c(pl_mAMY_r,file[2,])
		pl_lAMY_l=c(pl_lAMY_l,file[3,])
		pl_mAMY_l=c(pl_mAMY_l,file[4,])
	}
	# for each drug csv
        for (b in 1:length(Del_ps1FPs)){
                filename=paste0(paste0(childfp,'ses-01','/'),Del_ps1FPs[b])
                file=read.table(filename,header=FALSE,sep=',')
                ps_lAMY_r=c(ps_lAMY_r,file[1,])
                ps_mAMY_r=c(ps_mAMY_r,file[2,])
                ps_lAMY_l=c(ps_lAMY_l,file[3,])
		ps_mAMY_l=c(ps_mAMY_l,file[4,])
        }
}
# combine into dataframe
pl_outDF=data.frame(unlist(pl_lAMY_r),unlist(pl_mAMY_r),unlist(pl_lAMY_l),unlist(pl_mAMY_l))
ps_outDF=data.frame(unlist(ps_lAMY_r),unlist(ps_mAMY_r),unlist(ps_lAMY_l),unlist(ps_mAMY_l))

# now save em out
saveRDS(pl_outDF,'~/ptrd_pl_OutMagDF.rds')
saveRDS(ps_outDF,'~/ptrd_ps_OutMagDF.rds')






