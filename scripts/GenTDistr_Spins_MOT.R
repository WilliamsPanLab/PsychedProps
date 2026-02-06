# mirror load in from main script, but use it to generate t stats on spun distribution
library(nlme)
numSpins=2000
# initialize t distribution
tdistr_BUP=rep(0,numSpins)
tdistr_Mag=rep(0,numSpins)
tdistr_Seg=rep(0,numSpins)
tdistr_TA=rep(0,numSpins)

for (n in 1:numSpins){
	rs1fp=paste0('/scratch/users/apines/data/rs1_propsMerged_Spin_',n,'_MOT.csv')	
	rs2fp=paste0('/scratch/users/apines/data/rs2_propsMerged_Spin_',n,'_MOT.csv')	
	gamblingfp=paste0('/scratch/users/apines/data/gambling_propsMerged_Spin_',n,'_MOT.csv')	
	wmfp=paste0('/scratch/users/apines/data/wm_propsMerged_Spin_',n,'_MOT.csv')	
	# prop angles
	rs1=read.csv(rs1fp,header=F)
	rs2=read.csv(rs2fp,header=F)
	gambling=read.csv(gamblingfp,header=F)
	wm=read.csv(wmfp,header=F)	
	# set task
	rs1$Task='rs'
	rs2$Task='rs2'
	gambling$Task='gambling'
	wm$Task='wm'
	# subj 4 is out
	rs1=rs1[-c(4),]
	rs2=rs2[-c(4),]
	gambling=gambling[-c(4),]
	wm=wm[-c(4),]
	# manually pair columns as sep. observations of baseline, placebo, 80, 120mg
	rs1bv=data.frame(cbind(rs1$V1,rs1$V2,rs1$V3,rs1$V4,rs1$V17))
	rs1p=data.frame(cbind(rs1$V5,rs1$V6,rs1$V7,rs1$V8,rs1$V18))
	rs1m1=data.frame(cbind(rs1$V9,rs1$V10,rs1$V11,rs1$V12,rs1$V19))
	rs1m2=data.frame(cbind(rs1$V13,rs1$V14,rs1$V15,rs1$V16,rs1$V20))
	colnames(rs1bv)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs')
	colnames(rs1p)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs')
	colnames(rs1m1)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs')
	colnames(rs1m2)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs')

	rs2bv=data.frame(cbind(rs2$V1,rs2$V2,rs2$V3,rs2$V4,rs2$V17))
	rs2p=data.frame(cbind(rs2$V5,rs2$V6,rs2$V7,rs2$V8,rs2$V18))
	rs2m1=data.frame(cbind(rs2$V9,rs2$V10,rs2$V11,rs2$V12,rs2$V19))
	rs2m2=data.frame(cbind(rs2$V13,rs2$V14,rs2$V15,rs2$V16,rs2$V20))
	colnames(rs2bv)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs')
	colnames(rs2p)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs')
	colnames(rs2m1)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs')
	colnames(rs2m2)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs')

	gamblingbv=data.frame(cbind(gambling$V1,gambling$V2,gambling$V3,gambling$V4,gambling$V17))
	gamblingp=data.frame(cbind(gambling$V5,gambling$V6,gambling$V7,gambling$V8,gambling$V18))
	gamblingm1=data.frame(cbind(gambling$V9,gambling$V10,gambling$V11,gambling$V12,gambling$V19))
	gamblingm2=data.frame(cbind(gambling$V13,gambling$V14,gambling$V15,gambling$V16,gambling$V20))
	colnames(gamblingbv)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs')
	colnames(gamblingp)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs')
	colnames(gamblingm1)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs')
	colnames(gamblingm2)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs')

	wmbv=data.frame(cbind(wm$V1,wm$V2,wm$V3,wm$V4,wm$V17))
	wmp=data.frame(cbind(wm$V5,wm$V6,wm$V7,wm$V8,wm$V18))
	wmm1=data.frame(cbind(wm$V9,wm$V10,wm$V11,wm$V12,wm$V19))
	wmm2=data.frame(cbind(wm$V13,wm$V14,wm$V15,wm$V16,wm$V20))
	colnames(wmbv)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs')
	colnames(wmp)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs')
	colnames(wmm1)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs')
	colnames(wmm2)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs')
	# get subject IDs from this rds
	alff=readRDS('~/OutPlacDrug_alff.rds')
	rs1bv$Subjects=alff$SubjID
	rs1p$Subjects=alff$SubjID
	rs1m1$Subjects=alff$SubjID
	rs1m2$Subjects=alff$SubjID
	rs2bv$Subjects=alff$SubjID
	rs2p$Subjects=alff$SubjID
	rs2m1$Subjects=alff$SubjID
	rs2m2$Subjects=alff$SubjID
	gamblingbv$Subjects=alff$SubjID
	gamblingp$Subjects=alff$SubjID
	gamblingm1$Subjects=alff$SubjID
	gamblingm2$Subjects=alff$SubjID
	wmbv$Subjects=alff$SubjID
	wmp$Subjects=alff$SubjID
	wmm1$Subjects=alff$SubjID
	wmm2$Subjects=alff$SubjID
	# add in task (rs to be made equivalent after motion merge)
	rs1bv$Task='rs'
	rs1p$Task='rs'
	rs1m1$Task='rs'
	rs1m2$Task='rs'
	rs2bv$Task='rs2'
	rs2p$Task='rs2'
	rs2m1$Task='rs2'
	rs2m2$Task='rs2'
	gamblingbv$Task='gambling'
	gamblingp$Task='gambling'
	gamblingm1$Task='gambling'
	gamblingm2$Task='gambling'
	wmbv$Task='wm'
	wmp$Task='wm'
	wmm1$Task='wm'
	wmm2$Task='wm'
	# add in dosage
	rs1bv$Dosage='baseline'
	rs1p$Dosage='Placebo'
	rs1m1$Dosage='80mg'
	rs1m2$Dosage='120mg'
	rs2bv$Dosage='baseline'
	rs2p$Dosage='Placebo'
	rs2m1$Dosage='80mg'
	rs2m2$Dosage='120mg'
	gamblingbv$Dosage='baseline'
	gamblingp$Dosage='Placebo'
	gamblingm1$Dosage='80mg'
	gamblingm2$Dosage='120mg'
	wmbv$Dosage='baseline'
	wmp$Dosage='Placebo'
	wmm1$Dosage='80mg'
	wmm2$Dosage='120mg'
	# combine all
	allScans=rbind(rs1bv,rs1p,rs1m1,rs1m2,rs2bv,rs2p,rs2m1,rs2m2,gamblingbv,gamblingp,gamblingm1,gamblingm2,wmbv,wmp,wmm1,wmm2)
	# read in motion
	mot=read.csv('~/MDMA_spikes_summary.csv')
	# motion merge
	mergedDf=merge(mot,allScans,by=c('Subjects','Task','Dosage'))
	mergedDf$Drug=0
	mergedDf$Drug[mergedDf$Dosage=="120mg"]=1
	mergedDf$Drug[mergedDf$Dosage=="80mg"]=1
	mergedDf$Drug=as.factor(mergedDf$Drug)
	mergedDf$Subjects<-as.factor(mergedDf$Subjects)
	mergedDf$Dosage<-as.factor(mergedDf$Dosage)
	mergedDfProps=mergedDf
	# remove data that needs to be removed (subs 6 and 10, <250 TRs)
	mergedDf=mergedDf[mergedDf$Subjects!='sub-MDMA006',]
	mergedDf=mergedDf[mergedDf$Subjects!='sub-MDMA010',]
	mergedDf=mergedDf[mergedDf$RemTRs>250,]
	mergedDf$Task[mergedDf$Task=='rs2']='rs'
	mergedDf$Task=as.factor(mergedDf$Task)
	mergedDf_clean=mergedDf[mergedDf$Dosage!='baseline',]
	# fit null model
	fit_lme <- lme(TDProp1 ~ Drug + RemTRs + MeanFD+Task, random = ~ 1 | Subjects, data = mergedDf_clean)
	summaryLME<-summary(fit_lme)
	tdistr_BUP[n]=summaryLME$tTable['Drug1','t-value']
	
}

print('done with spun angles')

	###############
	#### MAGNITUDES

for (n in 1:numSpins){
        rs1fp=paste0('/scratch/users/apines/data/rs1_MagsMerged_Spin_',n,'_MOT.csv')
        rs2fp=paste0('/scratch/users/apines/data/rs2_MagsMerged_Spin_',n,'_MOT.csv')
        gamblingfp=paste0('/scratch/users/apines/data/gambling_MagsMerged_Spin_',n,'_MOT.csv')
        wmfp=paste0('/scratch/users/apines/data/wm_MagsMerged_Spin_',n,'_MOT.csv')
        # prop angles
        rs1=read.csv(rs1fp,header=F)
        rs2=read.csv(rs2fp,header=F)
        gambling=read.csv(gamblingfp,header=F)
        wm=read.csv(wmfp,header=F)
        # set task
        rs1$Task='rs'
        rs2$Task='rs2'
        gambling$Task='gambling'
        wm$Task='wm'
        # subj 4 is out
        rs1=rs1[-c(4),]
        rs2=rs2[-c(4),]
        gambling=gambling[-c(4),]
        wm=wm[-c(4),]
        # manually pair columns as sep. observations of baseline, placebo, 80, 120mg
        rs1bv=data.frame(cbind(rs1$V1,rs1$V2,rs1$V3,rs1$V4,rs1$V17))
        rs1p=data.frame(cbind(rs1$V5,rs1$V6,rs1$V7,rs1$V8,rs1$V18))
        rs1m1=data.frame(cbind(rs1$V9,rs1$V10,rs1$V11,rs1$V12,rs1$V19))
        rs1m2=data.frame(cbind(rs1$V13,rs1$V14,rs1$V15,rs1$V16,rs1$V20))
        colnames(rs1bv)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs')
        colnames(rs1p)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs')
        colnames(rs1m1)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs')
        colnames(rs1m2)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs')

        rs2bv=data.frame(cbind(rs2$V1,rs2$V2,rs2$V3,rs2$V4,rs2$V17))
        rs2p=data.frame(cbind(rs2$V5,rs2$V6,rs2$V7,rs2$V8,rs2$V18))
        rs2m1=data.frame(cbind(rs2$V9,rs2$V10,rs2$V11,rs2$V12,rs2$V19))
        rs2m2=data.frame(cbind(rs2$V13,rs2$V14,rs2$V15,rs2$V16,rs2$V20))
        colnames(rs2bv)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs')
        colnames(rs2p)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs')
        colnames(rs2m1)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs')
        colnames(rs2m2)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs')

        gamblingbv=data.frame(cbind(gambling$V1,gambling$V2,gambling$V3,gambling$V4,gambling$V17))
        gamblingp=data.frame(cbind(gambling$V5,gambling$V6,gambling$V7,gambling$V8,gambling$V18))
        gamblingm1=data.frame(cbind(gambling$V9,gambling$V10,gambling$V11,gambling$V12,gambling$V19))
        gamblingm2=data.frame(cbind(gambling$V13,gambling$V14,gambling$V15,gambling$V16,gambling$V20))
        colnames(gamblingbv)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs')
        colnames(gamblingp)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs')
        colnames(gamblingm1)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs')
        colnames(gamblingm2)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs')

        wmbv=data.frame(cbind(wm$V1,wm$V2,wm$V3,wm$V4,wm$V17))
        wmp=data.frame(cbind(wm$V5,wm$V6,wm$V7,wm$V8,wm$V18))
        wmm1=data.frame(cbind(wm$V9,wm$V10,wm$V11,wm$V12,wm$V19))
        wmm2=data.frame(cbind(wm$V13,wm$V14,wm$V15,wm$V16,wm$V20))
        colnames(wmbv)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs')
        colnames(wmp)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs')
        colnames(wmm1)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs')
        colnames(wmm2)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs')
        # get subject IDs from this rds
        alff=readRDS('~/OutPlacDrug_alff.rds')
        rs1bv$Subjects=alff$SubjID
        rs1p$Subjects=alff$SubjID
        rs1m1$Subjects=alff$SubjID
        rs1m2$Subjects=alff$SubjID
        rs2bv$Subjects=alff$SubjID
        rs2p$Subjects=alff$SubjID
        rs2m1$Subjects=alff$SubjID
        rs2m2$Subjects=alff$SubjID
        gamblingbv$Subjects=alff$SubjID
        gamblingp$Subjects=alff$SubjID
        gamblingm1$Subjects=alff$SubjID
        gamblingm2$Subjects=alff$SubjID
        wmbv$Subjects=alff$SubjID
        wmp$Subjects=alff$SubjID
        wmm1$Subjects=alff$SubjID
        wmm2$Subjects=alff$SubjID
        # add in task (rs to be made equivalent after motion merge)
        rs1bv$Task='rs'
        rs1p$Task='rs'
        rs1m1$Task='rs'
        rs1m2$Task='rs'
        rs2bv$Task='rs2'
        rs2p$Task='rs2'
        rs2m1$Task='rs2'
        rs2m2$Task='rs2'
        gamblingbv$Task='gambling'
        gamblingp$Task='gambling'
        gamblingm1$Task='gambling'
        gamblingm2$Task='gambling'
        wmbv$Task='wm'
        wmp$Task='wm'
        wmm1$Task='wm'
        wmm2$Task='wm'
        # add in dosage
        rs1bv$Dosage='baseline'
        rs1p$Dosage='Placebo'
        rs1m1$Dosage='80mg'
        rs1m2$Dosage='120mg'
        rs2bv$Dosage='baseline'
        rs2p$Dosage='Placebo'
        rs2m1$Dosage='80mg'
        rs2m2$Dosage='120mg'
        gamblingbv$Dosage='baseline'
        gamblingp$Dosage='Placebo'
        gamblingm1$Dosage='80mg'
        gamblingm2$Dosage='120mg'
        wmbv$Dosage='baseline'
        wmp$Dosage='Placebo'
        wmm1$Dosage='80mg'
        wmm2$Dosage='120mg'
        # combine all
        allScans=rbind(rs1bv,rs1p,rs1m1,rs1m2,rs2bv,rs2p,rs2m1,rs2m2,gamblingbv,gamblingp,gamblingm1,gamblingm2,wmbv,wmp,wmm1,wmm2)
        # read in motion
        mot=read.csv('~/MDMA_spikes_summary.csv')
        # motion merge
        mergedDf=merge(mot,allScans,by=c('Subjects','Task','Dosage'))
        mergedDf$Drug=0
        mergedDf$Drug[mergedDf$Dosage=="120mg"]=1
        mergedDf$Drug[mergedDf$Dosage=="80mg"]=1
        mergedDf$Drug=as.factor(mergedDf$Drug)
        mergedDf$Subjects<-as.factor(mergedDf$Subjects)
        mergedDf$Dosage<-as.factor(mergedDf$Dosage)
        mergedDfProps=mergedDf
        # remove data that needs to be removed (subs 6 and 10, <250 TRs)
        mergedDf=mergedDf[mergedDf$Subjects!='sub-MDMA006',]
        mergedDf=mergedDf[mergedDf$Subjects!='sub-MDMA010',]
        mergedDf=mergedDf[mergedDf$RemTRs>250,]
        mergedDf$Task[mergedDf$Task=='rs2']='rs'
        mergedDf$Task=as.factor(mergedDf$Task)
        mergedDf_clean=mergedDf[mergedDf$Dosage!='baseline',]
        # fit null model
        fit_lme <- lme(TDProp1 ~ Drug + RemTRs + MeanFD+Task, random = ~ 1 | Subjects, data = mergedDf_clean)
        summaryLME<-summary(fit_lme)
        tdistr_Mag[n]=summaryLME$tTable['Drug1','t-value']

}

print('done with DMN mag')

	############
	#### DMN SEG

for (n in 1:numSpins){
        rs1fp=paste0('/scratch/users/apines/data/rs1_MOTSegMerged_Spin_',n,'.csv')
        rs2fp=paste0('/scratch/users/apines/data/rs2_MOTSegMerged_Spin_',n,'.csv')
        gamblingfp=paste0('/scratch/users/apines/data/gambling_MOTSegMerged_Spin_',n,'.csv')
        wmfp=paste0('/scratch/users/apines/data/wm_MOTSegMerged_Spin_',n,'.csv')
        # prop angles
        rs1=read.csv(rs1fp,header=F)
        rs2=read.csv(rs2fp,header=F)
        gambling=read.csv(gamblingfp,header=F)
        wm=read.csv(wmfp,header=F)
        # set task
        rs1$Task='rs'
        rs2$Task='rs2'
        gambling$Task='gambling'
        wm$Task='wm'
        # subj 4 is out
        rs1=rs1[-c(4),]
        rs2=rs2[-c(4),]
        gambling=gambling[-c(4),]
        wm=wm[-c(4),]
        # manually pair columns as sep. observations of baseline, placebo, 80, 120mg
        rs1bv=data.frame(cbind(rs1$V1,rs1$V2,rs1$V3,rs1$V4,rs1$V17))
        rs1p=data.frame(cbind(rs1$V5,rs1$V6,rs1$V7,rs1$V8,rs1$V18))
        rs1m1=data.frame(cbind(rs1$V9,rs1$V10,rs1$V11,rs1$V12,rs1$V19))
        rs1m2=data.frame(cbind(rs1$V13,rs1$V14,rs1$V15,rs1$V16,rs1$V20))
        colnames(rs1bv)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs')
        colnames(rs1p)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs')
        colnames(rs1m1)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs')
        colnames(rs1m2)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs')

        rs2bv=data.frame(cbind(rs2$V1,rs2$V2,rs2$V3,rs2$V4,rs2$V17))
        rs2p=data.frame(cbind(rs2$V5,rs2$V6,rs2$V7,rs2$V8,rs2$V18))
        rs2m1=data.frame(cbind(rs2$V9,rs2$V10,rs2$V11,rs2$V12,rs2$V19))
        rs2m2=data.frame(cbind(rs2$V13,rs2$V14,rs2$V15,rs2$V16,rs2$V20))
        colnames(rs2bv)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs')
        colnames(rs2p)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs')
        colnames(rs2m1)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs')
        colnames(rs2m2)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs')

        gamblingbv=data.frame(cbind(gambling$V1,gambling$V2,gambling$V3,gambling$V4,gambling$V17))
        gamblingp=data.frame(cbind(gambling$V5,gambling$V6,gambling$V7,gambling$V8,gambling$V18))
        gamblingm1=data.frame(cbind(gambling$V9,gambling$V10,gambling$V11,gambling$V12,gambling$V19))
        gamblingm2=data.frame(cbind(gambling$V13,gambling$V14,gambling$V15,gambling$V16,gambling$V20))
        colnames(gamblingbv)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs')
        colnames(gamblingp)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs')
        colnames(gamblingm1)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs')
        colnames(gamblingm2)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs')

        wmbv=data.frame(cbind(wm$V1,wm$V2,wm$V3,wm$V4,wm$V17))
        wmp=data.frame(cbind(wm$V5,wm$V6,wm$V7,wm$V8,wm$V18))
        wmm1=data.frame(cbind(wm$V9,wm$V10,wm$V11,wm$V12,wm$V19))
        wmm2=data.frame(cbind(wm$V13,wm$V14,wm$V15,wm$V16,wm$V20))
        colnames(wmbv)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs')
        colnames(wmp)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs')
        colnames(wmm1)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs')
        colnames(wmm2)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs')
        # get subject IDs from this rds
        alff=readRDS('~/OutPlacDrug_alff.rds')
        rs1bv$Subjects=alff$SubjID
        rs1p$Subjects=alff$SubjID
        rs1m1$Subjects=alff$SubjID
        rs1m2$Subjects=alff$SubjID
        rs2bv$Subjects=alff$SubjID
        rs2p$Subjects=alff$SubjID
        rs2m1$Subjects=alff$SubjID
        rs2m2$Subjects=alff$SubjID
        gamblingbv$Subjects=alff$SubjID
        gamblingp$Subjects=alff$SubjID
        gamblingm1$Subjects=alff$SubjID
        gamblingm2$Subjects=alff$SubjID
        wmbv$Subjects=alff$SubjID
        wmp$Subjects=alff$SubjID
        wmm1$Subjects=alff$SubjID
        wmm2$Subjects=alff$SubjID
        # add in task (rs to be made equivalent after motion merge)
        rs1bv$Task='rs'
        rs1p$Task='rs'
        rs1m1$Task='rs'
        rs1m2$Task='rs'
        rs2bv$Task='rs2'
        rs2p$Task='rs2'
        rs2m1$Task='rs2'
        rs2m2$Task='rs2'
        gamblingbv$Task='gambling'
        gamblingp$Task='gambling'
        gamblingm1$Task='gambling'
        gamblingm2$Task='gambling'
        wmbv$Task='wm'
        wmp$Task='wm'
        wmm1$Task='wm'
        wmm2$Task='wm'
        # add in dosage
        rs1bv$Dosage='baseline'
        rs1p$Dosage='Placebo'
        rs1m1$Dosage='80mg'
        rs1m2$Dosage='120mg'
        rs2bv$Dosage='baseline'
        rs2p$Dosage='Placebo'
        rs2m1$Dosage='80mg'
        rs2m2$Dosage='120mg'
        gamblingbv$Dosage='baseline'
        gamblingp$Dosage='Placebo'
        gamblingm1$Dosage='80mg'
        gamblingm2$Dosage='120mg'
        wmbv$Dosage='baseline'
        wmp$Dosage='Placebo'
        wmm1$Dosage='80mg'
        wmm2$Dosage='120mg'
        # combine all
        allScans=rbind(rs1bv,rs1p,rs1m1,rs1m2,rs2bv,rs2p,rs2m1,rs2m2,gamblingbv,gamblingp,gamblingm1,gamblingm2,wmbv,wmp,wmm1,wmm2)
        # read in motion
        mot=read.csv('~/MDMA_spikes_summary.csv')
        # motion merge
        mergedDf=merge(mot,allScans,by=c('Subjects','Task','Dosage'))
        mergedDf$Drug=0
        mergedDf$Drug[mergedDf$Dosage=="120mg"]=1
        mergedDf$Drug[mergedDf$Dosage=="80mg"]=1
        mergedDf$Drug=as.factor(mergedDf$Drug)
        mergedDf$Subjects<-as.factor(mergedDf$Subjects)
        mergedDf$Dosage<-as.factor(mergedDf$Dosage)
        mergedDfProps=mergedDf
        # remove data that needs to be removed (subs 6 and 10, <250 TRs)
        mergedDf=mergedDf[mergedDf$Subjects!='sub-MDMA006',]
        mergedDf=mergedDf[mergedDf$Subjects!='sub-MDMA010',]
        mergedDf=mergedDf[mergedDf$RemTRs>250,]
        mergedDf$Task[mergedDf$Task=='rs2']='rs'
        mergedDf$Task=as.factor(mergedDf$Task)
        mergedDf_clean=mergedDf[mergedDf$Dosage!='baseline',]
        # fit null model
        fit_lme <- lme(TDProp1 ~ Drug + RemTRs + MeanFD+Task, random = ~ 1 | Subjects, data = mergedDf_clean)
        summaryLME<-summary(fit_lme)
        tdistr_Seg[n]=summaryLME$tTable['Drug1','t-value']

}

print('done with dmnseg')

	################
	#### DMN AUTOCOR

for (n in 1:numSpins){
        rs1fp=paste0('/scratch/users/apines/data/rs1_propsMerged_Spin_',n,'_MOT.csv')
        rs2fp=paste0('/scratch/users/apines/data/rs2_propsMerged_Spin_',n,'_MOT.csv')
        gamblingfp=paste0('/scratch/users/apines/data/gambling_propsMerged_Spin_',n,'_MOT.csv')
        wmfp=paste0('/scratch/users/apines/data/wm_propsMerged_Spin_',n,'_MOT.csv')
        # prop angles
        rs1=read.csv(rs1fp,header=F)
        rs2=read.csv(rs2fp,header=F)
        gambling=read.csv(gamblingfp,header=F)
        wm=read.csv(wmfp,header=F)
        # set task
        rs1$Task='rs'
        rs2$Task='rs2'
        gambling$Task='gambling'
        wm$Task='wm'
        # subj 4 is out
        rs1=rs1[-c(4),]
        rs2=rs2[-c(4),]
        gambling=gambling[-c(4),]
        wm=wm[-c(4),]
        # manually pair columns as sep. observations of baseline, placebo, 80, 120mg
        rs1bv=data.frame(cbind(rs1$V1,rs1$V2,rs1$V3,rs1$V4,rs1$V17))
        rs1p=data.frame(cbind(rs1$V5,rs1$V6,rs1$V7,rs1$V8,rs1$V18))
        rs1m1=data.frame(cbind(rs1$V9,rs1$V10,rs1$V11,rs1$V12,rs1$V19))
        rs1m2=data.frame(cbind(rs1$V13,rs1$V14,rs1$V15,rs1$V16,rs1$V20))
        colnames(rs1bv)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs')
        colnames(rs1p)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs')
        colnames(rs1m1)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs')
        colnames(rs1m2)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs')

        rs2bv=data.frame(cbind(rs2$V1,rs2$V2,rs2$V3,rs2$V4,rs2$V17))
        rs2p=data.frame(cbind(rs2$V5,rs2$V6,rs2$V7,rs2$V8,rs2$V18))
        rs2m1=data.frame(cbind(rs2$V9,rs2$V10,rs2$V11,rs2$V12,rs2$V19))
        rs2m2=data.frame(cbind(rs2$V13,rs2$V14,rs2$V15,rs2$V16,rs2$V20))
        colnames(rs2bv)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs')
        colnames(rs2p)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs')
        colnames(rs2m1)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs')
        colnames(rs2m2)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs')

        gamblingbv=data.frame(cbind(gambling$V1,gambling$V2,gambling$V3,gambling$V4,gambling$V17))
        gamblingp=data.frame(cbind(gambling$V5,gambling$V6,gambling$V7,gambling$V8,gambling$V18))
        gamblingm1=data.frame(cbind(gambling$V9,gambling$V10,gambling$V11,gambling$V12,gambling$V19))
        gamblingm2=data.frame(cbind(gambling$V13,gambling$V14,gambling$V15,gambling$V16,gambling$V20))
        colnames(gamblingbv)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs')
        colnames(gamblingp)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs')
        colnames(gamblingm1)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs')
        colnames(gamblingm2)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs')

        wmbv=data.frame(cbind(wm$V1,wm$V2,wm$V3,wm$V4,wm$V17))
        wmp=data.frame(cbind(wm$V5,wm$V6,wm$V7,wm$V8,wm$V18))
        wmm1=data.frame(cbind(wm$V9,wm$V10,wm$V11,wm$V12,wm$V19))
        wmm2=data.frame(cbind(wm$V13,wm$V14,wm$V15,wm$V16,wm$V20))
        colnames(wmbv)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs')
        colnames(wmp)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs')
        colnames(wmm1)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs')
        colnames(wmm2)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs')
        # get subject IDs from this rds
        alff=readRDS('~/OutPlacDrug_alff.rds')
        rs1bv$Subjects=alff$SubjID
        rs1p$Subjects=alff$SubjID
        rs1m1$Subjects=alff$SubjID
        rs1m2$Subjects=alff$SubjID
        rs2bv$Subjects=alff$SubjID
        rs2p$Subjects=alff$SubjID
        rs2m1$Subjects=alff$SubjID
        rs2m2$Subjects=alff$SubjID
        gamblingbv$Subjects=alff$SubjID
        gamblingp$Subjects=alff$SubjID
        gamblingm1$Subjects=alff$SubjID
        gamblingm2$Subjects=alff$SubjID
        wmbv$Subjects=alff$SubjID
        wmp$Subjects=alff$SubjID
        wmm1$Subjects=alff$SubjID
        wmm2$Subjects=alff$SubjID
        # add in task (rs to be made equivalent after motion merge)
        rs1bv$Task='rs'
        rs1p$Task='rs'
        rs1m1$Task='rs'
        rs1m2$Task='rs'
        rs2bv$Task='rs2'
        rs2p$Task='rs2'
        rs2m1$Task='rs2'
        rs2m2$Task='rs2'
        gamblingbv$Task='gambling'
        gamblingp$Task='gambling'
        gamblingm1$Task='gambling'
        gamblingm2$Task='gambling'
        wmbv$Task='wm'
        wmp$Task='wm'
        wmm1$Task='wm'
        wmm2$Task='wm'
        # add in dosage
        rs1bv$Dosage='baseline'
        rs1p$Dosage='Placebo'
        rs1m1$Dosage='80mg'
        rs1m2$Dosage='120mg'
        rs2bv$Dosage='baseline'
        rs2p$Dosage='Placebo'
        rs2m1$Dosage='80mg'
        rs2m2$Dosage='120mg'
        gamblingbv$Dosage='baseline'
        gamblingp$Dosage='Placebo'
        gamblingm1$Dosage='80mg'
        gamblingm2$Dosage='120mg'
        wmbv$Dosage='baseline'
        wmp$Dosage='Placebo'
        wmm1$Dosage='80mg'
        wmm2$Dosage='120mg'
        # combine all
        allScans=rbind(rs1bv,rs1p,rs1m1,rs1m2,rs2bv,rs2p,rs2m1,rs2m2,gamblingbv,gamblingp,gamblingm1,gamblingm2,wmbv,wmp,wmm1,wmm2)
        # read in motion
        mot=read.csv('~/MDMA_spikes_summary.csv')
        # motion merge
        mergedDf=merge(mot,allScans,by=c('Subjects','Task','Dosage'))
        mergedDf$Drug=0
        mergedDf$Drug[mergedDf$Dosage=="120mg"]=1
        mergedDf$Drug[mergedDf$Dosage=="80mg"]=1
        mergedDf$Drug=as.factor(mergedDf$Drug)
        mergedDf$Subjects<-as.factor(mergedDf$Subjects)
        mergedDf$Dosage<-as.factor(mergedDf$Dosage)
        mergedDfProps=mergedDf
        # remove data that needs to be removed (subs 6 and 10, <250 TRs)
        mergedDf=mergedDf[mergedDf$Subjects!='sub-MDMA006',]
        mergedDf=mergedDf[mergedDf$Subjects!='sub-MDMA010',]
        mergedDf=mergedDf[mergedDf$RemTRs>250,]
        mergedDf$Task[mergedDf$Task=='rs2']='rs'
        mergedDf$Task=as.factor(mergedDf$Task)
        mergedDf_clean=mergedDf[mergedDf$Dosage!='baseline',]
        # fit null model
        fit_lme <- lme(TDProp1 ~ Drug + RemTRs + MeanFD+Task, random = ~ 1 | Subjects, data = mergedDf_clean)
        summaryLME<-summary(fit_lme)
        tdistr_TA[n]=summaryLME$tTable['Drug1','t-value']

}

print ('done with TA aggregation')

# save out
outspins=data.frame(tdistr_BUP,tdistr_Mag,tdistr_Seg,tdistr_TA)
saveRDS(outspins,'~/SpunTDistributions_MDMA_MOT.rds')
