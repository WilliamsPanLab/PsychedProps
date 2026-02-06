# mirror load in from main script, but use it to generate t stats on spun distribution
library(nlme)
numSpins=2000
# initialize t distribution
tdistr_BUP=rep(0,numSpins)
tdistr_Mag=rep(0,numSpins)
tdistr_Seg=rep(0,numSpins)
tdistr_TA=rep(0,numSpins)

for (n in 1:numSpins){
	LSD=read.csv(paste0('/scratch/users/apines/data/lsd_propsMerged_Spin_',n,'_VIS.csv'))
	Thresh=200
	LSD=LSD[LSD$outDF_5>Thresh,]
	LSD=LSD[!is.na(LSD$outDF_3),]
	colnames(LSD)<-c('subj','sesh','PercBUP','task','RemainingTRs','meanFD')
	LSD$sesh <- relevel(as.factor(LSD$sesh), ref = "PCB")
	fit_lme <- lme(PercBUP ~ sesh + RemainingTRs + meanFD+task, random = ~ 1 | subj, data = LSD)
	summaryLME<-summary(fit_lme)
	tdistr_BUP[n]=summaryLME$tTable['seshLSD','t-value']
}
print('done with BUP')
for (n in 1:numSpins){
        LSD=read.csv(paste0('/scratch/users/apines/data/lsd_VISMagMerged_Spin_',n,'.csv'))
        LSD=LSD[LSD$outDF_5>Thresh,]
        LSD=LSD[!is.na(LSD$outDF_3),]
        colnames(LSD)<-c('subj','sesh','PercBUP','task','RemainingTRs','meanFD')
        LSD$sesh <- relevel(as.factor(LSD$sesh), ref = "PCB")
        fit_lme <- lme(PercBUP ~ sesh + RemainingTRs + meanFD+task, random = ~ 1 | subj, data = LSD)
        summaryLME<-summary(fit_lme)
        tdistr_Mag[n]=summaryLME$tTable['seshLSD','t-value']
}
print('done with mag')
for (n in 1:numSpins){
        LSD=read.csv(paste0('/scratch/users/apines/data/lsd_VISSegMerged_Spin_',n,'.csv'))
        LSD=LSD[LSD$outDF_5>Thresh,]
        LSD=LSD[!is.na(LSD$outDF_3),]
        colnames(LSD)<-c('subj','sesh','PercBUP','task','RemainingTRs','meanFD')
        LSD$sesh <- relevel(as.factor(LSD$sesh), ref = "PCB")
        fit_lme <- lme(PercBUP ~ sesh + RemainingTRs + meanFD+task, random = ~ 1 | subj, data = LSD)
        summaryLME<-summary(fit_lme)
        tdistr_Seg[n]=summaryLME$tTable['seshLSD','t-value']
}
print('done with seg')
for (n in 1:numSpins){
        LSD=read.csv(paste0('/scratch/users/apines/data/lsd_TAMerged_Spin_',n,'_VIS.csv'))
        LSD=LSD[LSD$outDF_5>Thresh,]
        LSD=LSD[!is.na(LSD$outDF_3),]
        colnames(LSD)<-c('subj','sesh','PercBUP','task','RemainingTRs','meanFD')
        LSD$sesh <- relevel(as.factor(LSD$sesh), ref = "PCB")
        fit_lme <- lme(PercBUP ~ sesh + RemainingTRs + meanFD+task, random = ~ 1 | subj, data = LSD)
        summaryLME<-summary(fit_lme)
        tdistr_TA[n]=summaryLME$tTable['seshLSD','t-value']
}
print('done with TA')
outspins=data.frame(tdistr_BUP,tdistr_Mag,tdistr_Seg,tdistr_TA)
saveRDS(outspins,'~/SpunTDistributions_LSD_VIS.rds')
