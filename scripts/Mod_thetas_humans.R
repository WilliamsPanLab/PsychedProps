################ libraries
library(nlme)
library(bpngreg)
################ load in MDMA
# baseline
rs1bv=read.csv('/oak/stanford/groups/leanew1/users/apines/data/rs1_thetas_L_bv.csv',header=F)
rs2bv=read.csv('/oak/stanford/groups/leanew1/users/apines/data/rs2_thetas_L_bv.csv',header=F)
emobv=read.csv('/oak/stanford/groups/leanew1/users/apines/data/emotion_thetas_L_bv.csv',header=F)
gamblingbv=read.csv('/oak/stanford/groups/leanew1/users/apines/data/gambling_thetas_L_bv.csv',header=F)
wmbv=read.csv('/oak/stanford/groups/leanew1/users/apines/data/wm_thetas_L_bv.csv',header=F)
rs1bv$Task='rs'
rs2bv$Task='rs'
emobv$Task='emotion'
gamblingbv$Task='gambling'
wmbv$Task='wm'
# remove subject 4
rs1bv=rs1bv[-c(4),]
rs2bv=rs2bv[-c(4),]
emobv=emobv[-c(4),]
gamblingbv=gamblingbv[-c(4),]
wmbv=wmbv[-c(4),]
# set last column of each to RemTRs
colnames(rs1bv)[ncol(rs1bv)]='RemTRs'
colnames(rs2bv)[ncol(rs2bv)]='RemTRs'
colnames(emobv)[ncol(emobv)]='RemTRs'
colnames(gamblingbv)[ncol(gamblingbv)]='RemTRs'
colnames(wmbv)[ncol(wmbv)]='RemTRs'
# get subject IDs and mot from this rds
subjInfo=readRDS('~/OutPlacDrug_alff.rds');
rs1bv$Subjects=subjInfo$SubjID
rs2bv$Subjects=subjInfo$SubjID
emobv$Subjects=subjInfo$SubjID
gamblingbv$Subjects=subjInfo$SubjID
wmbv$Subjects=subjInfo$SubjID
# Dosage as 'baseline'
rs1bv$Dosage='baseline'
rs2bv$Dosage='baseline'
emobv$Dosage='baseline'
gamblingbv$Dosage='baseline'
wmbv$Dosage='baseline'
# placebo
rs1pl=read.csv('/oak/stanford/groups/leanew1/users/apines/data/rs1_thetas_L_pcb.csv',header=F)
rs2pl=read.csv('/oak/stanford/groups/leanew1/users/apines/data/rs2_thetas_L_pcb.csv',header=F)
emopl=read.csv('/oak/stanford/groups/leanew1/users/apines/data/emotion_thetas_L_pcb.csv',header=F)
gamblingpl=read.csv('/oak/stanford/groups/leanew1/users/apines/data/gambling_thetas_L_pcb.csv',header=F)
wmpl=read.csv('/oak/stanford/groups/leanew1/users/apines/data/wm_thetas_L_pcb.csv',header=F)
rs1pl$Task='rs'
rs2pl$Task='rs'
emopl$Task='emotion'
gamblingpl$Task='gambling'
wmpl$Task='wm'
# remove subject 4
rs1pl=rs1pl[-c(4),]
rs2pl=rs2pl[-c(4),]
emopl=emopl[-c(4),]
gamblingpl=gamblingpl[-c(4),]
wmpl=wmpl[-c(4),]
# set last column of each to RemTRs
colnames(rs1pl)[ncol(rs1pl)]='RemTRs'
colnames(rs2pl)[ncol(rs2pl)]='RemTRs'
colnames(emopl)[ncol(emopl)]='RemTRs'
colnames(gamblingpl)[ncol(gamblingpl)]='RemTRs'
colnames(wmpl)[ncol(wmpl)]='RemTRs'
rs1pl$Subjects=subjInfo$SubjID
rs2pl$Subjects=subjInfo$SubjID
emopl$Subjects=subjInfo$SubjID
gamblingpl$Subjects=subjInfo$SubjID
wmpl$Subjects=subjInfo$SubjID
rs1pl$Dosage='Placebo'
rs2pl$Dosage='Placebo'
emopl$Dosage='Placebo'
gamblingpl$Dosage='Placebo'
wmpl$Dosage='Placebo'
# M1
rs1m1=read.csv('/oak/stanford/groups/leanew1/users/apines/data/rs1_thetas_L_m1.csv',header=F)
rs2m1=read.csv('/oak/stanford/groups/leanew1/users/apines/data/rs2_thetas_L_m1.csv',header=F)
emom1=read.csv('/oak/stanford/groups/leanew1/users/apines/data/emotion_thetas_L_m1.csv',header=F)
gamblingm1=read.csv('/oak/stanford/groups/leanew1/users/apines/data/gambling_thetas_L_m1.csv',header=F)
wmm1=read.csv('/oak/stanford/groups/leanew1/users/apines/data/wm_thetas_L_m1.csv',header=F)
rs1m1$Task='rs'
rs2m1$Task='rs'
emom1$Task='emotion'
gamblingm1$Task='gambling'
wmm1$Task='wm'
# remove subject 4
rs1m1=rs1m1[-c(4),]
rs2m1=rs2m1[-c(4),]
emom1=emom1[-c(4),]
gamblingm1=gamblingm1[-c(4),]
wmm1=wmm1[-c(4),]
# set last column of each to RemTRs
colnames(rs1m1)[ncol(rs1m1)]='RemTRs'
colnames(rs2m1)[ncol(rs2m1)]='RemTRs'
colnames(emom1)[ncol(emom1)]='RemTRs'
colnames(gamblingm1)[ncol(gamblingm1)]='RemTRs'
colnames(wmm1)[ncol(wmm1)]='RemTRs'
rs1m1$Subjects=subjInfo$SubjID
rs2m1$Subjects=subjInfo$SubjID
emom1$Subjects=subjInfo$SubjID
gamblingm1$Subjects=subjInfo$SubjID
wmm1$Subjects=subjInfo$SubjID
rs1m1$Dosage='80mg'
rs2m1$Dosage='80mg'
emom1$Dosage='80mg'
gamblingm1$Dosage='80mg'
wmm1$Dosage='80mg'
# M2
rs1m2=read.csv('/oak/stanford/groups/leanew1/users/apines/data/rs1_thetas_L_m2.csv',header=F)
rs2m2=read.csv('/oak/stanford/groups/leanew1/users/apines/data/rs2_thetas_L_m2.csv',header=F)
emom2=read.csv('/oak/stanford/groups/leanew1/users/apines/data/emotion_thetas_L_m2.csv',header=F)
gamblingm2=read.csv('/oak/stanford/groups/leanew1/users/apines/data/gambling_thetas_L_m2.csv',header=F)
wmm2=read.csv('/oak/stanford/groups/leanew1/users/apines/data/wm_thetas_L_m2.csv',header=F)
rs1m2$Task='rs'
rs2m2$Task='rs'
emom2$Task='emotion'
gamblingm2$Task='gambling'
wmm2$Task='wm'
# remove subject 4
rs1m2=rs1m2[-c(4),]
rs2m2=rs2m2[-c(4),]
emom2=emom2[-c(4),]
gamblingm2=gamblingm2[-c(4),]
wmm2=wmm2[-c(4),]
# set last column of each to RemTRs
colnames(rs1m2)[ncol(rs1m2)]='RemTRs'
colnames(rs2m2)[ncol(rs2m2)]='RemTRs'
colnames(emom2)[ncol(emom2)]='RemTRs'
colnames(gamblingm2)[ncol(gamblingm2)]='RemTRs'
colnames(wmm2)[ncol(wmm2)]='RemTRs'
rs1m2$Subjects=subjInfo$SubjID
rs2m2$Subjects=subjInfo$SubjID
emom2$Subjects=subjInfo$SubjID
gamblingm2$Subjects=subjInfo$SubjID
wmm2$Subjects=subjInfo$SubjID
rs1m2$Dosage='120mg'
rs2m2$Dosage='120mg'
emom2$Dosage='120mg'
gamblingm2$Dosage='120mg'
wmm2$Dosage='120mg'
# combine all MDMA data
allScans_m=rbind(rs1bv, rs2bv, emobv, gamblingbv, wmbv,
		rs1pl, rs2pl, emopl, gamblingpl, wmpl,
		rs1m1, rs2m1, emom1, gamblingm1, wmm1,
		rs1m2, rs2m2, emom2, gamblingm2, wmm2)

# merge in motion
mot=read.csv('~/MDMA_spikes_summary.csv')
mergedMDMA=merge(mot, allScans_m, by=c('Subjects', 'Task', 'Dosage'))
mergedDfMDMA$Drug=0
mergedDfMDMA$Drug[mergedDf$Dosage=="120mg"]=1
mergedDfMDMA$Drug[mergedDf$Dosage=="80mg"]=1
mergedDfMDMA$Drug=as.factor(mergedDf$Drug)
mergedDfMDMA$Subjects=as.factor(mergedDf$Subjects)

#################### load in LSD
# placebo
rs1pl=read.csv('/oak/stanford/groups/leanew1/users/apines/data/lsd_thetas_rs1_PCB_L.csv',header=F)
rs2pl=read.csv('/oak/stanford/groups/leanew1/users/apines/data/lsd_thetas_rs2_PCB_L.csv',header=F)
muspl=read.csv('/oak/stanford/groups/leanew1/users/apines/data/lsd_thetas_mus_PCB_L.csv',header=F)
# LSD
rs1lsd=read.csv('/oak/stanford/groups/leanew1/users/apines/data/lsd_thetas_rs1_LSD_L.csv',header=F)
rs2lsd=read.csv('/oak/stanford/groups/leanew1/users/apines/data/lsd_thetas_rs2_LSD_L.csv',header=F)
muslsd=read.csv('/oak/stanford/groups/leanew1/users/apines/data/lsd_thetas_mus_LSD_L.csv',header=F)

#################### load in psil

#################### model each face

#################### save out stats
