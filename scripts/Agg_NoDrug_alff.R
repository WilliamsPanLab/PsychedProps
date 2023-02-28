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

# initialize df out here with subj names, each ROI name x2 for drug nondrug
df<-data.frame(subjList)
colnames(df)<-'SubjID'

# for each subj
for (s in 1:length(subjvec)){
	# initialize these vectors sep. for each subject for averaging
	pl_lAMY_r=NULL
	pl_mAMY_r=NULL
	pl_lAMY_l=NULL
	pl_mAMY_l=NULL
	md_lAMY_r=NULL
	md_mAMY_r=NULL
	md_lAMY_l=NULL
	md_mAMY_l=NULL
	# and for thalamus
	pl_THA-VPm-rh=NULL
	pl_THA-VPl-rh=NULL
	pl_THA-VAi-rh=NULL
	pl_THA-VAs-rh=NULL
	pl_THA-DAm-rh=NULL
	pl_THA-DAl-rh=NULL
	pl_THA-VPm-lh=NULL
	pl_THA-VPl-lh=NULL
	pl_THA-VAi-lh=NULL
	pl_THA-VAs-lh=NULL
	pl_THA-DAm-lh=NULL
	pl_THA-DAl-lh=NULL
	md_THA-VPm-rh=NULL
	md_THA-VPl-rh=NULL
	md_THA-VAi-rh=NULL
	md_THA-VAs-rh=NULL
	md_THA-DAm-rh=NULL
	md_THA-DAl-rh=NULL
	md_THA-VPm-lh=NULL
	md_THA-VPl-lh=NULL
	md_THA-VAi-lh=NULL
	md_THA-VAs-lh=NULL
	md_THA-DAm-lh=NULL
	md_THA-DAl-lh=NULL
	# get subj ID
	subj=subjList[s]
	# plop subject ID (confirmatory)
	df$SubjID[s]=subj
	# get output filepath
	childfp=paste0('/oak/stanford/groups/leanew1/users/apines/data/p50/',subj,'/') 
	Del_baselineFPs=list.files(path=(paste0(childfp,p1[s],'/')),pattern="_alff.txt")
	Del_pFPs=list.files(path=(paste0(childfp,p2[s],'/')),pattern="_alff.txt")
	Del_m1FPs=list.files(path=(paste0(childfp,m1[s],'/')),pattern="_alff.txt")
	Del_m2FPs=list.files(path=(paste0(childfp,m2[s],'/')),pattern="_alff.txt")
	# load in each non-drug csv
	for (b in 1:length(Del_baselineFPs)){
		filename=paste0(paste0(childfp,p1[s],'/'),Del_baselineFPs[b])
		file=read.table(filename)
		pl_lAMY_r=c(pl_lAMY_r,file[19,])
		pl_mAMY_r=c(pl_mAMY_r,file[20,])
		pl_lAMY_l=c(pl_lAMY_l,file[44,])
		pl_mAMY_l=c(pl_mAMY_l,file[45,])
		pl_THA-VPm-rh=c(pl_THA-VPm-rh,file[5,])
		pl_THA-VPl-rh=c(pl_THA-VPl-rh,file[6,])
		pl_THA-VAi-rh=c(pl_THA-VAi-rh,file[7,])
		pl_THA-VAs-rh=c(pl_THA-VAs-rh,file[8,])
		pl_THA-DAm-rh=c(pl_THA-DAm-rh,file[9,])
		pl_THA-DAl-rh=c(pl_THA-DAl-rh,file[10,])
		pl_THA-VPm-lh=c(pl_THA-VPm-lh,file[30,])
		pl_THA-VPl-lh=c(pl_THA-VPl-lh,file[31,])
		pl_THA-VAi-lh=c(pl_THA-VAi-lh,file[32,])
		pl_THA-VAs-lh=c(pl_THA-VAs-lh,file[33,])
		pl_THA-DAm-lh=c(pl_THA-DAm-lh,file[34,])
		pl_THA-DAl-lh=c(pl_THA-DAl-lh,file[35,])
	}
	for (b in 1:length(Del_pFPs)){
		filename=paste0(paste0(childfp,p2[s],'/'),Del_pFPs[b])
		file=read.table(filename)
                pl_lAMY_r=c(pl_lAMY_r,file[19,])
                pl_mAMY_r=c(pl_mAMY_r,file[20,])
                pl_lAMY_l=c(pl_lAMY_l,file[44,])
                pl_mAMY_l=c(pl_mAMY_l,file[45,])
        }
	# for each drug csv
        for (b in 1:length(Del_m1FPs)){
                filename=paste0(paste0(childfp,m1[s],'/'),Del_m1FPs[b])
                file=read.table(filename)
                md_lAMY_r=c(md_lAMY_r,file[19,])
                md_mAMY_r=c(md_mAMY_r,file[20,])
                md_lAMY_l=c(md_lAMY_l,file[44,])
		md_mAMY_l=c(md_mAMY_l,file[45,])
        }
        for (b in 1:length(Del_m2FPs)){
                filename=paste0(paste0(childfp,m2[s],'/'),Del_m2FPs[b])
                file=read.table(filename)
                md_lAMY_r=c(md_lAMY_r,file[19,])
                md_mAMY_r=c(md_mAMY_r,file[20,])
                md_lAMY_l=c(md_lAMY_l,file[44,])
                md_mAMY_l=c(md_mAMY_l,file[45,])
        }
}
# combine into dataframe
outDF=data.frame(unlist(pl_lAMY_r),unlist(pl_mAMY_r),unlist(pl_lAMY_l),unlist(pl_mAMY_l),unlist(md_lAMY_r),unlist(md_mAMY_r),unlist(md_lAMY_l),unlist(md_mAMY_l))

# now save em out
saveRDS(outDF,'~/OutDelay_alff.rds')






