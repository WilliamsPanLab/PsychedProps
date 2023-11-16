library(nlme)
# take iteration number in as only argument (1-6 for left hemisphere)
iteration = as.numeric(commandArgs(trailingOnly=TRUE))

# Check if the argument is provided
if (is.na(iteration) || iteration < 1 || iteration > 7) {
	stop("I'm intended for iterating over the faces in chunks: 1, 2, 3, 4, 5, or 6 7. Tell me which chunk you want or find yourself another script.")
}

# Calculate the range of iterations based on the input
start_iteration <- (iteration - 1) * 271 + 1
end_iteration <- iteration * 271

# aggregated covariate information
covs=readRDS('/oak/stanford/groups/leanew1/users/apines/data/P50_cleaned_df.rds')
# change rs to rs1 to harmonize with later script
covs$Task[covs$Task=='rs']='rs1'
# extract range of vertices to be covered in this run from VertBin
Lfaces=1800
Rfaces=1897

### MANUALLY RECREATE SUBJ LIST
subjPrefix=rep('sub-MDMA0',17)
subjSuffix=c('01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17')
subjList=paste0(subjPrefix,subjSuffix)

# initialize face-level vectors: keeping 0s in the slots untouched by this run to verify allocation later
Tg_R=rep(0,Rfaces)
pg_R=rep(0,Rfaces)
Tw_R=rep(0,Rfaces)
pw_R=rep(0,Rfaces)

# create dataframe for output
df<-data.frame(subjList)
colnames(df)<-'SubjID'

# initialize iterable face value column - to avoid storing this data for multiple faces over this upcoming loop
df$DMNProp_pl=rep(NA,length(subjList))
df$DMNProp_m1=rep(NA,length(subjList))
df$DMNProp_m2=rep(NA,length(subjList))
df$DMNProp_pl_FD=rep(NA,length(subjList))
df$DMNProp_m1_FD=rep(NA,length(subjList))
df$DMNProp_m2_FD=rep(NA,length(subjList))
df$DMNProp_pl_RemTRs=rep(NA,length(subjList))
df$DMNProp_m1_RemTRs=rep(NA,length(subjList))
df$DMNProp_m2_RemTRs=rep(NA,length(subjList))
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

# saveout sub-session-dose correspondence
write.table(data.frame(p1,p2,m1,m2),'~/subjSeshDoseCorresp.csv',quote=F,col.names=F)

# subjvec: no 4
subjvec=c(1,2,3,5,7,8,9,11,12,13,14,15,16,17)
tasks=c('rs1','rs2','wm','gambling')
# for each Right hemi face
for (f in start_iteration:end_iteration){
	print(f)
	ctmasterdf=NULL
	# for each task
	for (task in tasks){
		CovRowsTask=covs[covs$Task==task,]
		# load in DMN props iteratively for each subj
		for (s in subjvec){
			subj=subjList[s]
			# get covariate rows corresponding to this subject
			CovRowsSubj=CovRowsTask[CovRowsTask$Subjects==subj,]
			#file paths
			Lp2FP=paste0('/scratch/users/apines/data/mdma/',subj,'/',p2[s],'/',subj,'_',p2[s],'_',task,'_Prop_TS_dmn_R.csv')
			Lm1FP=paste0('/scratch/users/apines/data/mdma/',subj,'/',m1[s],'/',subj,'_',m1[s],'_',task,'_Prop_TS_dmn_R.csv')
			Lm2FP=paste0('/scratch/users/apines/data/mdma/',subj,'/',m2[s],'/',subj,'_',m2[s],'_',task,'_Prop_TS_dmn_R.csv')
			
			# if Lp2FP exists
			if (file.exists(Lp2FP)){
				# load in dat data
				Lp2=read.csv(Lp2FP,header=F)
				# extract DMN prop from placebo
				Prop_pl=mean(unlist(Lp2[f,]))
				df$DMNProp_pl[s]=Prop_pl
				# extract FD
				df$DMNProp_pl_FD[s]=CovRowsSubj$MeanFD[CovRowsSubj$Dosage=='Placebo']
				# extract remaining TRs
				df$DMNProp_pl_RemTRs[s]=CovRowsSubj$RemTRs[CovRowsSubj$Dosage=='Placebo']
			}
			
			# if Lm1FP exists
			if (file.exists(Lm1FP)){
				# load in dat data
				Lm1=read.csv(Lm1FP,header=F)
				# extract DMN prop from placebo
				Prop_m1=mean(unlist(Lm1[f,]))
				df$DMNProp_m1[s]=Prop_m1
				# extract FD
				df$DMNProp_m1_FD[s]=CovRowsSubj$MeanFD[CovRowsSubj$Dosage=='80mg']
				# extract remaining TRs
				df$DMNProp_m1_RemTRs[s]=CovRowsSubj$RemTRs[CovRowsSubj$Dosage=='80mg']
			}

			# if Lm2FP exists
			if (file.exists(Lm2FP)){
				# load in dat data
				Lm2=read.csv(Lm2FP,header=F)
				# extract DMN prop from placebo
				Prop_m2=mean(unlist(Lm2[f,]))
				df$DMNProp_m2[s]=Prop_m2
				# extract FD
				df$DMNProp_m2_FD[s]=CovRowsSubj$MeanFD[CovRowsSubj$Dosage=='120mg']
				# extract remaining TRs
				df$DMNProp_m2_RemTRs[s]=CovRowsSubj$RemTRs[CovRowsSubj$Dosage=='120mg']
			}
		}
		# make a seperate df for pl, m1, and m2
		pldf=data.frame(df$SubjID,df$DMNProp_pl,df$DMNProp_pl_FD,df$DMNProp_pl_RemTRs)
		colnames(pldf)=c('SubjID','DMNProp','FD','RemTRs')
		pldf$Drug=0
		m1df=data.frame(df$SubjID,df$DMNProp_m1,df$DMNProp_m1_FD,df$DMNProp_m1_RemTRs)
		colnames(m1df)=c('SubjID','DMNProp','FD','RemTRs')
		m1df$Drug=1
		m2df=data.frame(df$SubjID,df$DMNProp_m2,df$DMNProp_m2_FD,df$DMNProp_m2_RemTRs)
		colnames(m2df)=c('SubjID','DMNProp','FD','RemTRs')
		m2df$Drug=1
		# convert to long-format df
		masterdf=rbind(pldf,m1df,m2df)
		# apply same TR threshold as main analyses
		masterdf=masterdf[masterdf$RemTRs>250,]
		# remove NA values (true missing)
		masterdf=na.omit(masterdf)
		# add task variable
		masterdf$Task=task
		# add it to cross-task master df
		ctmasterdf=rbind(ctmasterdf,masterdf)
		}
	# code rs as one task
        ctmasterdf$Task[ctmasterdf$Task=='rs1']='rs'
        ctmasterdf$Task[ctmasterdf$Task=='rs2']='rs'
        # set rs to reference
        ctmasterdf$Task<-as.factor(ctmasterdf$Task)
        ctmasterdf$Task <- relevel(ctmasterdf$Task, ref = "rs")
        # mixed model with subject intercept for Prop and covarying for motion and remtrs in testing Drug condition
        mmodel=lme(DMNProp ~ Drug + Task+ FD + RemTRs, random = ~ 1 | SubjID, data = ctmasterdf)
        # extract t stat
        Tg_R[f]=summary(mmodel)$tTable[3,4]
        pg_R[f]=summary(mmodel)$tTable[3,5]
        Tw_R[f]=summary(mmodel)$tTable[4,4]
        pw_R[f]=summary(mmodel)$tTable[4,5]
}
# write out t stats and ps
write.csv(Tg_R,paste0('/oak/stanford/groups/leanew1/users/apines/results/DMNProps_Tstats_R_Taskg_',as.character(iteration),'.csv'))
write.csv(pg_R,paste0('/oak/stanford/groups/leanew1/users/apines/results/DMNProps_ps_R_Taskg_',as.character(iteration),'.csv'))
write.csv(Tw_R,paste0('/oak/stanford/groups/leanew1/users/apines/results/DMNProps_Tstats_R_Taskw_',as.character(iteration),'.csv'))
write.csv(pw_R,paste0('/oak/stanford/groups/leanew1/users/apines/results/DMNProps_ps_R_Taskw_',as.character(iteration),'.csv'))
