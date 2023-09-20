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
# save these for looking up motion later
Num_p1=rep(0,length(subjList))
Num_p2=p2
Num_m1=m1
Num_m2=m2

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

# initialize row counter (first row is initialization values
RowCount=2

# initialize masterdf
drug<-'Init'
Pulses<-0
Subject<-'Init'
scan<-'Init'
FD<-0.00
ROI<-'roi'

# make ROI vectors: we are going to loop over them instead of populating them sep.
ROIList=c('pl_lAMY_r','pl_mAMY_r','pl_lAMY_l','pl_mAMY_l','HIP-head-m-rh','HIP-head-l-rh','HIP-body-rh','HIP-tail-rh','HIP-head-m-lh','HIP-head-l-lh','HIP-body-lh','HIP-tail-lh','17Networks_LH_SalVentAttnB_PFCmp_1','17Networks_RH_SalVentAttnB_PFCmp_1','17Networks_LH_DefaultB_PFCv_1','17Networks_LH_SalVentAttnA_Ins_2','17Networks_LH_SalVentAttnA_Ins_1','17Networks_RH_DefaultB_PFCv_1','17Networks_RH_SalVentAttnA_Ins_1')
ROINumb=seq(1,19)

# for each subj
for (s in 1:length(subjvec)){
	subj=subjList[s]
	# get output filepath
	childfp=paste0('/oak/stanford/groups/leanew1/users/apines/data/p50/',subj,'/') 
	Del_baselineFPs=list.files(path=(paste0(childfp,p1[s],'/')),pattern="_PulseCount.csv")
	Del_pFPs=list.files(path=(paste0(childfp,p2[s],'/')),pattern="_PulseCount.csv")
	Del_m1FPs=list.files(path=(paste0(childfp,m1[s],'/')),pattern="_PulseCount.csv")
	Del_m2FPs=list.files(path=(paste0(childfp,m2[s],'/')),pattern="_PulseCount.csv")
	# load in each non-drug csv
	for (b in 1:length(Del_baselineFPs)){
		# get filename
		filename=paste0(paste0(childfp,p1[s],'/'),Del_baselineFPs[b])
		file=read.table(filename,header=FALSE,sep=',')
		# for each ROI
		for (r in 1:19){
			# populate masterdf vectors
			drug[RowCount]<-'No'
			Pulses[RowCount]<-file[1,ROINumb[r]]
			Subject[RowCount]<-subj
			# determine scan type to populate scan type 
			task=strsplit(filename,'_')[[1]][2]
			scan[RowCount]<-task
			# FD is eveyr entry in row 2 just take 1
			FD[RowCount]<-file[2,2]
			ROI[RowCount]=ROIList[r]
			# new row for next ROI
			RowCount=RowCount+1
	
		}
	}
	# placebo tp 1
	for (b in 1:length(Del_pFPs)){
		# get filename
		filename=paste0(paste0(childfp,p2[s],'/'),Del_pFPs[b])
		file=read.table(filename,header=FALSE,sep=',')
		for (r in 1:19){
			# populate masterdf vectors
			drug[RowCount]<-'No'
			Pulses[RowCount]<-file[1,ROINumb[r]]
			Subject[RowCount]<-subj
			# determine scan type to populate scan type 
			task=strsplit(filename,'_')[[1]][2]
			scan[RowCount]<-task
			# FD is eveyr entry in row 2 just take 1
			FD[RowCount]<-file[2,2]
			ROI[RowCount]=ROIList[r]
			# new row for next ROI
			RowCount=RowCount+1
		}
        }
	# md tp 1
	# for each drug csv
        for (b in 1:length(Del_m1FPs)){
        	# get filename
                filename=paste0(paste0(childfp,m1[s],'/'),Del_m1FPs[b])
		file=read.table(filename,header=FALSE,sep=',')
		for (r in 1:19){
			# populate masterdf vectors
			drug[RowCount]<-'MDMA'
			Pulses[RowCount]<-file[1,ROINumb[r]]
			Subject[RowCount]<-subj
			# determine scan type to populate scan type 
			task=strsplit(filename,'_')[[1]][2]
			scan[RowCount]<-task
			# FD is eveyr entry in row 2 just take 1
			FD[RowCount]<-file[2,2]
			ROI[RowCount]=ROIList[r]
			# new row for next ROI
			RowCount=RowCount+1
		}
        }
	# md tp2
        for (b in 1:length(Del_m2FPs)){
        	# get filename
                filename=paste0(paste0(childfp,m2[s],'/'),Del_m2FPs[b])
		file=read.table(filename,header=FALSE,sep=',')	
		for (r in 1:19){
			# populate masterdf vectors
			drug[RowCount]<-'MDMA'
			Pulses[RowCount]<-file[1,ROINumb[r]]
			Subject[RowCount]<-subj
			# determine scan type to populate scan type 
			task=strsplit(filename,'_')[[1]][2]
			scan[RowCount]<-task
			# FD is eveyr entry in row 2 just take 1
			FD[RowCount]<-file[2,2]
			ROI[RowCount]=ROIList[r]
			# new row for next ROI
			RowCount=RowCount+1
		}
        }
}
# combine into masterdf
masterdf<-data.frame(cbind(drug,Pulses,Subject,scan,FD,ROI))

# now save em out
saveRDS(masterdf,'~/OutPulseDF.rds')






