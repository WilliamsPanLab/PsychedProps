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
PowerInfo<-c(0,0.00938967,0.01877934,0.02816901,0.03755869,0.04694836,0.05633803,0.0657277,0.07511737,0.08450704,0.09389671,0.10328638,0.11267606,0.12206573,0.1314554,0.14084507,0.15023474,0.15962441,0.16901408,0.17840376,0.18779343,0.1971831,0.20657277,0.21596244,0.22535211,0.23474178,0.24413146,0.25352113,0.2629108,0.27230047,0.28169014,0.29107981,0.30046948,0.30985915,0.31924883,0.3286385,0.33802817,0.34741784,0.35680751,0.36619718,0.37558685,0.38497653,0.3943662,0.40375587,0.41314554,0.42253521,0.43192488,0.44131455,0.45070423,0.4600939,0.46948357,0.47887324,0.48826291,0.49765258,0.50704225,0.51643192,0.5258216,0.53521127,0.54460094,0.55399061,0.56338028,0.57276995,0.58215962,0.5915493,0.60093897,0.61032864,0.61971831,0.62910798,0.63849765,0.64788732,0.657277,0.66666667,0.67605634,0.68544601,0.69483568,0.70422535)
Subject<-'Init'
scan<-'Init'
FD<-0.00
ROI<-'roi'

# initialize big output matrix
PowerDf=data.frame(matrix(nrow=1,ncol=81))

colnames(PowerDf)[1:5]=c('Drug','Subject','scan','ROI','FD')
PowerDf[1,6:81]=PowerInfo

# make ROI vectors: we are going to loop over them instead of populating them sep.
ROIList=c('pl_lAMY_r','pl_mAMY_r','pl_lAMY_l','pl_mAMY_l','HIP-head-m-rh','HIP-head-l-rh','HIP-body-rh','HIP-tail-rh','HIP-head-m-lh','HIP-head-l-lh','HIP-body-lh','HIP-tail-lh','17Networks_LH_SalVentAttnB_PFCmp_1','17Networks_RH_SalVentAttnB_PFCmp_1','17Networks_LH_DefaultB_PFCv_1','17Networks_LH_SalVentAttnA_Ins_2','17Networks_LH_SalVentAttnA_Ins_1','17Networks_RH_DefaultB_PFCv_1','17Networks_RH_SalVentAttnA_Ins_1')
ROINumb=seq(1,19)

# for each subj
for (s in 1:length(subjvec)){
	subj=subjList[s]
	# get output filepath
	childfp=paste0('/oak/stanford/groups/leanew1/users/apines/data/p50/',subj,'/') 
	Del_baselineFPs=list.files(path=(paste0(childfp,p1[s],'/')),pattern="_PSDs.csv")
	Del_pFPs=list.files(path=(paste0(childfp,p2[s],'/')),pattern="_PSDs.csv")
	Del_m1FPs=list.files(path=(paste0(childfp,m1[s],'/')),pattern="_PSDs.csv")
	Del_m2FPs=list.files(path=(paste0(childfp,m2[s],'/')),pattern="_PSDs.csv")
	# load in each non-drug csv
	for (b in 1:length(Del_baselineFPs)){
		# get filename
		filename=paste0(paste0(childfp,p1[s],'/'),Del_baselineFPs[b])
		file=read.table(filename,header=FALSE,sep=',')
		# determine scan type to populate scan type 
		task=strsplit(filename,'_')[[1]][2]
		# conf file
		confFilepath=paste0('/oak/stanford/groups/leanew1/SHARED_DATASETS/private/p50/bids/data/derivatives/fmriprep-20.2.3/fmriprep/',subj,'/',p1[s],'/func/',subj,'_',p1[s],'_task-',task,'_acq-mb_dir-pe0_run-0_desc-confounds_timeseries.tsv')
		if (task=='rs2'){
			confFilepath=paste0('/oak/stanford/groups/leanew1/SHARED_DATASETS/private/p50/bids/data/derivatives/fmriprep-20.2.3/fmriprep/',subj,'/',p1[s],'/func/',subj,'_',p1[s],'_task-rs_acq-mb_dir-pe1_run-0_desc-confounds_timeseries.tsv')
		}
		# read in conf file to get fd
		confFile=read.delim(confFilepath,sep="\t")
		# calculate FD
		FD=mean(na.omit(as.numeric(confFile$framewise_displacement)))
		# for each ROI
		for (r in 1:19){
			# convert to vector and aggregate for this ROI
			PowerDf[RowCount,]=c('No',subj,task,ROIList[r],FD,file[,r])
			# new row for next ROI
			RowCount=RowCount+1
		}
	}
	# placebo tp 1
	for (b in 1:length(Del_pFPs)){
		# get filename
		filename=paste0(paste0(childfp,p2[s],'/'),Del_pFPs[b])
		file=read.table(filename,header=FALSE,sep=',')
		# determine scan type to populate scan type 
                task=strsplit(filename,'_')[[1]][2]
                # conf file
                confFilepath=paste0('/oak/stanford/groups/leanew1/SHARED_DATASETS/private/p50/bids/data/derivatives/fmriprep-20.2.3/fmriprep/',subj,'/',p2[s],'/func/',subj,'_',p2[s],'_task-',task,'_acq-mb_dir-pe0_run-0_desc-confounds_timeseries.tsv')
                if (task=='rs2'){
  			confFilepath=paste0('/oak/stanford/groups/leanew1/SHARED_DATASETS/private/p50/bids/data/derivatives/fmriprep-20.2.3/fmriprep/',subj,'/',p2[s],'/func/',subj,'_',p2[s],'_task-rs_acq-mb_dir-pe1_run-0_desc-confounds_timeseries.tsv')
		}
                # read in conf file to get fd
                confFile=read.delim(confFilepath,sep="\t")
                # calculate FD
		FD=mean(na.omit(as.numeric(confFile$framewise_displacement)))
		for (r in 1:19){
			# convert to vector and aggregate for this ROI
                        PowerDf[RowCount,]=c('No',subj,task,ROIList[r],FD,file[,r])
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
		# determine scan type to populate scan type 
		task=strsplit(filename,'_')[[1]][2]
		# conf file
		confFilepath=paste0('/oak/stanford/groups/leanew1/SHARED_DATASETS/private/p50/bids/data/derivatives/fmriprep-20.2.3/fmriprep/',subj,'/',m1[s],'/func/',subj,'_',m1[s],'_task-',task,'_acq-mb_dir-pe0_run-0_desc-confounds_timeseries.tsv')
		if (task=='rs2'){      
			        confFilepath=paste0('/oak/stanford/groups/leanew1/SHARED_DATASETS/private/p50/bids/data/derivatives/fmriprep-20.2.3/fmriprep/',subj,'/',m1[s],'/func/',subj,'_',m1[s],'_task-rs_acq-mb_dir-pe1_run-0_desc-confounds_timeseries.tsv')
		}       
		# read in conf file to get fd
		confFile=read.delim(confFilepath,sep="\t")
		# calculate FD
		FD=mean(na.omit(as.numeric(confFile$framewise_displacement)))
		for (r in 1:19){
			# convert to vector and aggregate for this ROI
			PowerDf[RowCount,]=c('MDMA',subj,task,ROIList[r],FD,file[,r])
			RowCount=RowCount+1
		}
        }
	# md tp2
        for (b in 1:length(Del_m2FPs)){
        	# get filename
                filename=paste0(paste0(childfp,m2[s],'/'),Del_m2FPs[b])
		file=read.table(filename,header=FALSE,sep=',')	
		# determine scan type to populate scan type 
		task=strsplit(filename,'_')[[1]][2]
		# conf file
		confFilepath=paste0('/oak/stanford/groups/leanew1/SHARED_DATASETS/private/p50/bids/data/derivatives/fmriprep-20.2.3/fmriprep/',subj,'/',m2[s],'/func/',subj,'_',m2[s],'_task-',task,'_acq-mb_dir-pe0_run-0_desc-confounds_timeseries.tsv')
		if (task=='rs2'){      
			        confFilepath=paste0('/oak/stanford/groups/leanew1/SHARED_DATASETS/private/p50/bids/data/derivatives/fmriprep-20.2.3/fmriprep/',subj,'/',m2[s],'/func/',subj,'_',m2[s],'_task-rs_acq-mb_dir-pe1_run-0_desc-confounds_timeseries.tsv')
		}       
		# read in conf file to get fd
		confFile=read.delim(confFilepath,sep="\t")
		# calculate FD
		FD=mean(na.omit(as.numeric(confFile$framewise_displacement)))
		for (r in 1:19){
			# convert to vector and aggregate for this ROI
			PowerDf[RowCount,]=c('MDMA',subj,task,ROIList[r],FD,file[,r])
			RowCount=RowCount+1
		}
        }
}
# combine into masterdf
masterdf<-data.frame(PowerDf)

# now save em out
saveRDS(masterdf,'~/OutPowerDF.rds')






