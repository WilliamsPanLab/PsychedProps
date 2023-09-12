### MANUALLY RECREATE SUBJ LIST
subjPrefix=rep('sub-MDMA0',15)
# no subj 4 or 15
subjSuffix=c('01','02','03','05','06','07','08','09','10','11','12','13','14','15','16','17')
subjList=paste0(subjPrefix,subjSuffix)

# load session_dose information
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
subjvec=c(1,2,3,5,6,7,8,9,10,11,12,13,14,15,16,17)

# initialize df out here with subj names, each ROI name x2 for drug nondrug
df<-data.frame(subjList)
colnames(df)<-'SubjID'
# and ROI_specific column names
df$pl_lAMY_r=NULL
df$pl_mAMY_r=NULL
df$pl_lAMY_l=NULL
df$pl_mAMY_l=NULL
df$md_lAMY_r=NULL
df$md_mAMY_r=NULL
df$md_lAMY_l=NULL
df$md_mAMY_l=NULL
df$pl_THA_VPm_rh=NULL
df$pl_THA_VPl_rh=NULL
df$pl_THA_VAi_rh=NULL
df$pl_THA_VAs_rh=NULL
df$pl_THA_DAm_rh=NULL
df$pl_THA_DAl_rh=NULL
df$pl_THA_VPm_lh=NULL
df$pl_THA_VPl_lh=NULL
df$pl_THA_VAi_lh=NULL
df$pl_THA_VAs_lh=NULL
df$pl_THA_DAm_lh=NULL
df$pl_THA_DAl_lh=NULL
df$md_THA_VPm_rh=NULL
df$md_THA_VPl_rh=NULL
df$md_THA_VAi_rh=NULL
df$md_THA_VAs_rh=NULL
df$md_THA_DAm_rh=NULL
df$md_THA_DAl_rh=NULL
df$md_THA_VPm_lh=NULL
df$md_THA_VPl_lh=NULL
df$md_THA_VAi_lh=NULL
df$md_THA_VAs_lh=NULL
df$md_THA_DAm_lh=NULL
df$md_THA_DAl_lh=NULL
df$md80_THA_DAm=0
df$md120_THA_DAm=0
df$bv_THA_DAm=0
df$pla_THA_DAm=0
# for each subj
for (s in 1:length(subjvec)){
	print(s)
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
	pl_THA_VPm_rh=NULL
	pl_THA_VPl_rh=NULL
	pl_THA_VAi_rh=NULL
	pl_THA_VAs_rh=NULL
	pl_THA_DAm_rh=NULL
	pl_THA_DAl_rh=NULL
	pl_THA_VPm_lh=NULL
	pl_THA_VPl_lh=NULL
	pl_THA_VAi_lh=NULL
	pl_THA_VAs_lh=NULL
	pl_THA_DAm_lh=NULL
	pl_THA_DAl_lh=NULL
	md_THA_VPm_rh=NULL
	md_THA_VPl_rh=NULL
	md_THA_VAi_rh=NULL
	md_THA_VAs_rh=NULL
	md_THA_DAm_rh=NULL
	md_THA_DAl_rh=NULL
	md_THA_VPm_lh=NULL
	md_THA_VPl_lh=NULL
	md_THA_VAi_lh=NULL
	md_THA_VAs_lh=NULL
	md_THA_DAm_lh=NULL
	md_THA_DAl_lh=NULL
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
	# load in each non_drug csv
	for (b in 1:2){
		filename=paste0(paste0(childfp,p1[s],'/'),Del_baselineFPs[b])
		file=read.table(filename)
		pl_lAMY_r=c(pl_lAMY_r,file[19,])
		pl_mAMY_r=c(pl_mAMY_r,file[20,])
		pl_lAMY_l=c(pl_lAMY_l,file[44,])
		pl_mAMY_l=c(pl_mAMY_l,file[45,])
		pl_THA_VPm_rh=c(pl_THA_VPm_rh,file[5,])
		pl_THA_VPl_rh=c(pl_THA_VPl_rh,file[6,])
		pl_THA_VAi_rh=c(pl_THA_VAi_rh,file[7,])
		pl_THA_VAs_rh=c(pl_THA_VAs_rh,file[8,])
		pl_THA_DAm_rh=c(pl_THA_DAm_rh,file[9,])
		pl_THA_DAl_rh=c(pl_THA_DAl_rh,file[10,])
		pl_THA_VPm_lh=c(pl_THA_VPm_lh,file[30,])
		pl_THA_VPl_lh=c(pl_THA_VPl_lh,file[31,])
		pl_THA_VAi_lh=c(pl_THA_VAi_lh,file[32,])
		pl_THA_VAs_lh=c(pl_THA_VAs_lh,file[33,])
		pl_THA_DAm_lh=c(pl_THA_DAm_lh,file[34,])
		pl_THA_DAl_lh=c(pl_THA_DAl_lh,file[35,])
		df$bv_THA_DAm[s]=df$bv_THA_DAm[s]+(file[9,]+file[34,])
	}
	for (b in 1:2){
		filename=paste0(paste0(childfp,p2[s],'/'),Del_pFPs[b])
		file=read.table(filename)
                pl_lAMY_r=c(pl_lAMY_r,file[19,])
                pl_mAMY_r=c(pl_mAMY_r,file[20,])
                pl_lAMY_l=c(pl_lAMY_l,file[44,])
                pl_mAMY_l=c(pl_mAMY_l,file[45,])
        	pl_lAMY_r=c(pl_lAMY_r,file[19,])
		pl_mAMY_r=c(pl_mAMY_r,file[20,])
		pl_lAMY_l=c(pl_lAMY_l,file[44,])
		pl_mAMY_l=c(pl_mAMY_l,file[45,])
		pl_THA_VPm_rh=c(pl_THA_VPm_rh,file[5,])
		pl_THA_VPl_rh=c(pl_THA_VPl_rh,file[6,])
		pl_THA_VAi_rh=c(pl_THA_VAi_rh,file[7,])
		pl_THA_VAs_rh=c(pl_THA_VAs_rh,file[8,])
		pl_THA_DAm_rh=c(pl_THA_DAm_rh,file[9,])
		pl_THA_DAl_rh=c(pl_THA_DAl_rh,file[10,])
		pl_THA_VPm_lh=c(pl_THA_VPm_lh,file[30,])
		pl_THA_VPl_lh=c(pl_THA_VPl_lh,file[31,])
		pl_THA_VAi_lh=c(pl_THA_VAi_lh,file[32,])
		pl_THA_VAs_lh=c(pl_THA_VAs_lh,file[33,])
		pl_THA_DAm_lh=c(pl_THA_DAm_lh,file[34,])
		pl_THA_DAl_lh=c(pl_THA_DAl_lh,file[35,])
		df$pla_THA_DAm[s]=df$pla_THA_DAm[s]+(file[9,]+file[34,])
	}
	# for each drug csv
        for (b in 1:2){
                filename=paste0(paste0(childfp,m1[s],'/'),Del_m1FPs[b])
                file=read.table(filename)
                md_lAMY_r=c(md_lAMY_r,file[19,])
                md_mAMY_r=c(md_mAMY_r,file[20,])
                md_lAMY_l=c(md_lAMY_l,file[44,])
		md_mAMY_l=c(md_mAMY_l,file[45,])
		md_lAMY_r=c(pl_lAMY_r,file[19,])
		md_mAMY_r=c(pl_mAMY_r,file[20,])
		md_lAMY_l=c(pl_lAMY_l,file[44,])
		md_mAMY_l=c(pl_mAMY_l,file[45,])
		md_THA_VPm_rh=c(md_THA_VPm_rh,file[5,])
		md_THA_VPl_rh=c(md_THA_VPl_rh,file[6,])
		md_THA_VAi_rh=c(md_THA_VAi_rh,file[7,])
		md_THA_VAs_rh=c(md_THA_VAs_rh,file[8,])
		md_THA_DAm_rh=c(md_THA_DAm_rh,file[9,])
		md_THA_DAl_rh=c(md_THA_DAl_rh,file[10,])
		md_THA_VPm_lh=c(md_THA_VPm_lh,file[30,])
		md_THA_VPl_lh=c(md_THA_VPl_lh,file[31,])
		md_THA_VAi_lh=c(md_THA_VAi_lh,file[32,])
		md_THA_VAs_lh=c(md_THA_VAs_lh,file[33,])
		md_THA_DAm_lh=c(md_THA_DAm_lh,file[34,])
		md_THA_DAl_lh=c(md_THA_DAl_lh,file[35,])
        	df$md80_THA_DAm[s]=df$md80_THA_DAm[s]+(file[9,]+file[34,])
	}
        for (b in 1:2){
                filename=paste0(paste0(childfp,m2[s],'/'),Del_m2FPs[b])
                file=read.table(filename)
                md_lAMY_r=c(md_lAMY_r,file[19,])
                md_mAMY_r=c(md_mAMY_r,file[20,])
                md_lAMY_l=c(md_lAMY_l,file[44,])
                md_mAMY_l=c(md_mAMY_l,file[45,])
		md_lAMY_r=c(pl_lAMY_r,file[19,])
		md_mAMY_r=c(pl_mAMY_r,file[20,])
		md_lAMY_l=c(pl_lAMY_l,file[44,])
		md_mAMY_l=c(pl_mAMY_l,file[45,])
		md_THA_VPm_rh=c(md_THA_VPm_rh,file[5,])
		md_THA_VPl_rh=c(md_THA_VPl_rh,file[6,])
		md_THA_VAi_rh=c(md_THA_VAi_rh,file[7,])
		md_THA_VAs_rh=c(md_THA_VAs_rh,file[8,])
		md_THA_DAm_rh=c(md_THA_DAm_rh,file[9,])
		md_THA_DAl_rh=c(md_THA_DAl_rh,file[10,])
		md_THA_VPm_lh=c(md_THA_VPm_lh,file[30,])
		md_THA_VPl_lh=c(md_THA_VPl_lh,file[31,])
		md_THA_VAi_lh=c(md_THA_VAi_lh,file[32,])
		md_THA_VAs_lh=c(md_THA_VAs_lh,file[33,])
		md_THA_DAm_lh=c(md_THA_DAm_lh,file[34,])
		md_THA_DAl_lh=c(md_THA_DAl_lh,file[35,])
		df$md120_THA_DAm[s]=df$md120_THA_DAm[s]+(file[9,]+file[34,])
	}
	# average each vector for each subject for placebo and drug sep.
	df$pl_lAMY_r[s]=mean(na.omit(pl_lAMY_r))
	df$pl_mAMY_r[s]=mean(na.omit(pl_mAMY_r))
	df$pl_lAMY_l[s]=mean(na.omit(pl_lAMY_l))
	df$pl_mAMY_l[s]=mean(na.omit(pl_mAMY_l))
	df$md_lAMY_r[s]=mean(na.omit(md_lAMY_r))
	df$md_mAMY_r[s]=mean(na.omit(md_mAMY_r))
	df$md_lAMY_l[s]=mean(na.omit(md_lAMY_l))
	df$md_mAMY_l[s]=mean(na.omit(md_mAMY_l))
	df$pl_THA_VPm_rh[s]=mean(na.omit(pl_THA_VPm_rh))
	df$pl_THA_VPl_rh[s]=mean(na.omit(pl_THA_VPl_rh))
	df$pl_THA_VAi_rh[s]=mean(na.omit(pl_THA_VAi_rh))
	df$pl_THA_VAs_rh[s]=mean(na.omit(pl_THA_VAs_rh))
	df$pl_THA_DAm_rh[s]=mean(na.omit(pl_THA_DAm_rh))
	df$pl_THA_DAl_rh[s]=mean(na.omit(pl_THA_DAl_rh))
	df$pl_THA_VPm_lh[s]=mean(na.omit(pl_THA_VPm_lh))
	df$pl_THA_VPl_lh[s]=mean(na.omit(pl_THA_VPl_lh))
	df$pl_THA_VAi_lh[s]=mean(na.omit(pl_THA_VAi_lh))
	df$pl_THA_VAs_lh[s]=mean(na.omit(pl_THA_VAs_lh))
	df$pl_THA_DAm_lh[s]=mean(na.omit(pl_THA_DAm_lh))
	df$pl_THA_DAl_lh[s]=mean(na.omit(pl_THA_DAl_lh))
	df$md_THA_VPm_rh[s]=mean(na.omit(md_THA_VPm_rh))
	df$md_THA_VPl_rh[s]=mean(na.omit(md_THA_VPl_rh))
	df$md_THA_VAi_rh[s]=mean(na.omit(md_THA_VAi_rh))
	df$md_THA_VAs_rh[s]=mean(na.omit(md_THA_VAs_rh))
	df$md_THA_DAm_rh[s]=mean(na.omit(md_THA_DAm_rh))
	df$md_THA_DAl_rh[s]=mean(na.omit(md_THA_DAl_rh))
	df$md_THA_VPm_lh[s]=mean(na.omit(md_THA_VPm_lh))
	df$md_THA_VPl_lh[s]=mean(na.omit(md_THA_VPl_lh))
	df$md_THA_VAi_lh[s]=mean(na.omit(md_THA_VAi_lh))
	df$md_THA_VAs_lh[s]=mean(na.omit(md_THA_VAs_lh))
	df$md_THA_DAm_lh[s]=mean(na.omit(md_THA_DAm_lh))
	df$md_THA_DAl_lh[s]=mean(na.omit(md_THA_DAl_lh))
}
# now save em out
saveRDS(df,'~/OutPlacDrug_alff.rds')






