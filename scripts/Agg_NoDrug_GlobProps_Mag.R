### MANUALLY RECREATE SUBJ LIST
subjPrefix=rep('sub-MDMA0',15)
# no subj 4 or 15
subjSuffix=c('01','02','03','05','06','07','08','09','10','11','12','13','14','16','17')
subjList=paste0(subjPrefix,subjSuffix)

# load session-dose information
seshDose=read.csv('~/MDMA_Dose_info_unblinded_first17.csv',stringsAsFactors=FALSE)
### make session indices for each subject
# baseline as placebo
p1=rep('ses-00',length(subjList))
# placebo
p2=NULL
# dose one
m1=NULL
# dose two
m2=NULL
# for each subj
for (s in 1:length(subjList)){
	doseInfo=seshDose[s,2:5]
	p2[s]=grep('Placebo',doseInfo)
	m1[s]=grep('80mg',doseInfo)
	m2[s]=grep('120mg',doseInfo)
}

# replace numeric indices with ses code
p2<-gsub(2,'ses-01',p2)
p2<-gsub(3,'ses-02',p2)
p2<-gsub(4,'ses-03',p2)
m1<-gsub(2,'ses-01',m1)
m1<-gsub(3,'ses-02',m1)
m1<-gsub(4,'ses-03',m1)
m2<-gsub(2,'ses-01',m2)
m2<-gsub(3,'ses-02',m2)
m2<-gsub(4,'ses-03',m2)

# subjvec: no 4 (no 15 either)
subjvec=c(1,2,3,5,6,7,8,9,10,11,12,13,14,16,17)

# create output vectors for group-level inference on delay in each ROI
pl_lAMY_r=NULL
pl_mAMY_r=NULL
pl_lAMY_l=NULL
pl_mAMY_l=NULL
md_lAMY_r=NULL
md_mAMY_r=NULL
md_lAMY_l=NULL
md_mAMY_l=NULL

# for each subj
for (s in 1:length(subjvec)){
	subj=subjList[s]
	# get output filepath
	childfp=paste0('/oak/stanford/groups/leanew1/users/apines/data/p50/',subj,'/') 
	Del_baselineFPs=list.files(path=(paste0(childfp,p1[s],'/')),pattern="MagMat.csv")
	Del_pFPs=list.files(path=(paste0(childfp,p2[s],'/')),pattern="MagMat.csv")
	Del_m1FPs=list.files(path=(paste0(childfp,m1[s],'/')),pattern="MagMat.csv")
	Del_m2FPs=list.files(path=(paste0(childfp,m2[s],'/')),pattern="MagMat.csv")
	# load in each non-drug csv
	for (b in 1:length(Del_baselineFPs)){
		filename=paste0(paste0(childfp,p1[s],'/'),Del_baselineFPs[b])
		file=read.table(filename,header=FALSE,sep=',')
		pl_lAMY_r=c(pl_lAMY_r,file[1,])
		pl_mAMY_r=c(pl_mAMY_r,file[2,])
		pl_lAMY_l=c(pl_lAMY_l,file[3,])
		pl_mAMY_l=c(pl_mAMY_l,file[4,])
	}
	for (b in 1:length(Del_pFPs)){
		filename=paste0(paste0(childfp,p2[s],'/'),Del_pFPs[b])
		file=read.table(filename,header=FALSE,sep=',')
                pl_lAMY_r=c(pl_lAMY_r,file[1,])
                pl_mAMY_r=c(pl_mAMY_r,file[2,])
                pl_lAMY_l=c(pl_lAMY_l,file[3,])
                pl_mAMY_l=c(pl_mAMY_l,file[4,])
        }
	# for each drug csv
        for (b in 1:length(Del_m1FPs)){
                filename=paste0(paste0(childfp,m1[s],'/'),Del_m1FPs[b])
                file=read.table(filename,header=FALSE,sep=',')
                md_lAMY_r=c(md_lAMY_r,file[1,])
                md_mAMY_r=c(md_mAMY_r,file[2,])
                md_lAMY_l=c(md_lAMY_l,file[3,])
		md_mAMY_l=c(md_mAMY_l,file[4,])
        }
        for (b in 1:length(Del_m2FPs)){
                filename=paste0(paste0(childfp,m2[s],'/'),Del_m2FPs[b])
                file=read.table(filename,header=FALSE,sep=',')
                md_lAMY_r=c(md_lAMY_r,file[1,])
                md_mAMY_r=c(md_mAMY_r,file[2,])
                md_lAMY_l=c(md_lAMY_l,file[3,])
                md_mAMY_l=c(md_mAMY_l,file[4,])
        }
}
# combine into dataframe
pl_outDF=data.frame(unlist(pl_lAMY_r),unlist(pl_mAMY_r),unlist(pl_lAMY_l),unlist(pl_mAMY_l))
md_outDF=data.frame(unlist(md_lAMY_r),unlist(md_mAMY_r),unlist(md_lAMY_l),unlist(md_mAMY_l))

# now save em out
saveRDS(pl_outDF,'~/pl_OutMagDF.rds')
saveRDS(md_outDF,'~/md_OutMagDF.rds')






