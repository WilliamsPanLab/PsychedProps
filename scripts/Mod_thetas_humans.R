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
# set task
rs1bv$Task='rs'
rs2bv$Task='rs'
emobv$Task='emotion'
gamblingbv$Task='gambling'
wmbv$Task='wm'
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
# set task
rs1pl$Task='rs'
rs2pl$Task='rs'
emopl$Task='emotion'
gamblingpl$Task='gambling'
wmpl$Task='wm'
# set subj ID
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
# set task
rs1m1$Task='rs'
rs2m1$Task='rs'
emom1$Task='emotion'
gamblingm1$Task='gambling'
wmm1$Task='wm'
# set subj id
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
# set task
rs1m2$Task='rs'
rs2m2$Task='rs'
emom2$Task='emotion'
gamblingm2$Task='gambling'
wmm2$Task='wm'
# set subj id
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
mergedMDMA$Drug=0
mergedMDMA$Drug[mergedMDMA$Dosage=="120mg"]=1
mergedMDMA$Drug[mergedMDMA$Dosage=="80mg"]=1
mergedMDMA$Drug=as.factor(mergedMDMA$Drug)
mergedMDMA$Subjects=as.factor(mergedMDMA$Subjects)
# QC threshold

#################### load in LSD
# placebo
rs1pl=read.csv('/oak/stanford/groups/leanew1/users/apines/data/lsd_thetas_rs1_PCB_L.csv',header=F)
rs2pl=read.csv('/oak/stanford/groups/leanew1/users/apines/data/lsd_thetas_rs2_PCB_L.csv',header=F)
muspl=read.csv('/oak/stanford/groups/leanew1/users/apines/data/lsd_thetas_mus_PCB_L.csv',header=F)
# LSD
rs1lsd=read.csv('/oak/stanford/groups/leanew1/users/apines/data/lsd_thetas_rs1_LSD_L.csv',header=F)
rs2lsd=read.csv('/oak/stanford/groups/leanew1/users/apines/data/lsd_thetas_rs2_LSD_L.csv',header=F)
muslsd=read.csv('/oak/stanford/groups/leanew1/users/apines/data/lsd_thetas_mus_LSD_L.csv',header=F)
# set meanFD 
colnames(rs1pl)[ncol(rs1pl)]='MeanFD'
colnames(rs2pl)[ncol(rs2pl)]='MeanFD'
colnames(muspl)[ncol(muspl)]='MeanFD'
colnames(rs1lsd)[ncol(rs1lsd)]='MeanFD'
colnames(rs2lsd)[ncol(rs2lsd)]='MeanFD'
colnames(muslsd)[ncol(muslsd)]='MeanFD'
# set remTRs
colnames(rs1pl)[ncol(rs1pl)-1]='RemTRs'
colnames(rs2pl)[ncol(rs2pl)-1]='RemTRs'
colnames(muspl)[ncol(muspl)-1]='RemTRs'
colnames(rs1lsd)[ncol(rs1lsd)-1]='RemTRs'
colnames(rs2lsd)[ncol(rs2lsd)-1]='RemTRs'
colnames(muslsd)[ncol(muslsd)-1]='RemTRs'
# set tasks
rs1pl$Task='rs'
rs2pl$Task='rs'
rs1lsd$Task='rs'
rs2lsd$Task='rs'
muspl$Task='mus'
muslsd$Task='mus'
# set Drug
rs1pl$Drug=0
rs2pl$Drug=0
rs1lsd$Drug=1
rs2lsd$Drug=1
muspl$Drug=0
muslsd$Drug=1
# set subjects: nice and easy in this dataset
rs1pl$Subjects=paste0('sub-LSD',seq(1:20))
rs2pl$Subjects=paste0('sub-LSD',seq(1:20))
muspl$Subjects=paste0('sub-LSD',seq(1:20))
rs1lsd$Subjects=paste0('sub-LSD',seq(1:20))
rs2lsd$Subjects=paste0('sub-LSD',seq(1:20))
muslsd$Subjects=paste0('sub-LSD',seq(1:20))
# merge them together
allScans_l=rbind(rs1pl, rs2pl, muspl,
                rs1lsd, rs2lsd, muslsd)
# QC threshold

# as.factor factors
allScans_l$Subjects=as.factor(allScans_l$Subjects)
allScans_l$Drug=as.factor(allScans_l$Drug)
#################### load in psil
rs1BV=read.csv('/oak/stanford/groups/leanew1/users/apines/data/rs1_Psil_thetas_L_bv.csv',header=F)
rs2BV=read.csv('/oak/stanford/groups/leanew1/users/apines/data/rs2_Psil_thetas_L_bv.csv',header=F)
rs3BV=read.csv('/oak/stanford/groups/leanew1/users/apines/data/rs3_Psil_thetas_L_bv.csv',header=F)
rs4BV=read.csv('/oak/stanford/groups/leanew1/users/apines/data/rs4_Psil_thetas_L_bv.csv',header=F)
rs5BV=read.csv('/oak/stanford/groups/leanew1/users/apines/data/rs5_Psil_thetas_L_bv.csv',header=F)
rs6BV=read.csv('/oak/stanford/groups/leanew1/users/apines/data/rs6_Psil_thetas_L_bv.csv',header=F)
# between drugs
rs1BW=read.csv('/oak/stanford/groups/leanew1/users/apines/data/rs1_Psil_thetas_L_p.csv',header=F)
rs2BW=read.csv('/oak/stanford/groups/leanew1/users/apines/data/rs2_Psil_thetas_L_p.csv',header=F)
rs3BW=read.csv('/oak/stanford/groups/leanew1/users/apines/data/rs3_Psil_thetas_L_p.csv',header=F)
rs4BW=read.csv('/oak/stanford/groups/leanew1/users/apines/data/rs4_Psil_thetas_L_p.csv',header=F)
rs5BW=read.csv('/oak/stanford/groups/leanew1/users/apines/data/rs5_Psil_thetas_L_p.csv',header=F)
rs6BW=read.csv('/oak/stanford/groups/leanew1/users/apines/data/rs6_Psil_thetas_L_p.csv',header=F)
# after drug
rs1AF=read.csv('/oak/stanford/groups/leanew1/users/apines/data/rs1_Psil_thetas_L_m1.csv',header=F)
rs2AF=read.csv('/oak/stanford/groups/leanew1/users/apines/data/rs2_Psil_thetas_L_m1.csv',header=F)
rs3AF=read.csv('/oak/stanford/groups/leanew1/users/apines/data/rs3_Psil_thetas_L_m1.csv',header=F)
rs4AF=read.csv('/oak/stanford/groups/leanew1/users/apines/data/rs4_Psil_thetas_L_m1.csv',header=F)
rs5AF=read.csv('/oak/stanford/groups/leanew1/users/apines/data/rs5_Psil_thetas_L_m1.csv',header=F)
rs6AF=read.csv('/oak/stanford/groups/leanew1/users/apines/data/rs6_Psil_thetas_L_m1.csv',header=F)
# during drug
rs1D=read.csv('/oak/stanford/groups/leanew1/users/apines/data/rs1_Psil_thetas_L_m2.csv',header=F)
rs2D=read.csv('/oak/stanford/groups/leanew1/users/apines/data/rs2_Psil_thetas_L_m2.csv',header=F)
rs3D=read.csv('/oak/stanford/groups/leanew1/users/apines/data/rs3_Psil_thetas_L_m2.csv',header=F)
rs4D=read.csv('/oak/stanford/groups/leanew1/users/apines/data/rs4_Psil_thetas_L_m2.csv',header=F)
rs5D=read.csv('/oak/stanford/groups/leanew1/users/apines/data/rs5_Psil_thetas_L_m2.csv',header=F)
rs6D=read.csv('/oak/stanford/groups/leanew1/users/apines/data/rs6_Psil_thetas_L_m2.csv',header=F)
# set last column of each to meanFD
colnames(rs1BV)[ncol(rs1BV)]='MeanFD'
colnames(rs2BV)[ncol(rs2BV)]='MeanFD'
colnames(rs3BV)[ncol(rs3BV)]='MeanFD'
colnames(rs4BV)[ncol(rs4BV)]='MeanFD'
colnames(rs5BV)[ncol(rs5BV)]='MeanFD'
colnames(rs6BV)[ncol(rs6BV)]='MeanFD'
colnames(rs1BW)[ncol(rs1BW)]='MeanFD'
colnames(rs2BW)[ncol(rs2BW)]='MeanFD'
colnames(rs3BW)[ncol(rs3BW)]='MeanFD'
colnames(rs4BW)[ncol(rs4BW)]='MeanFD'
colnames(rs5BW)[ncol(rs5BW)]='MeanFD'
colnames(rs6BW)[ncol(rs6BW)]='MeanFD'
colnames(rs1AF)[ncol(rs1AF)]='MeanFD'
colnames(rs2AF)[ncol(rs2AF)]='MeanFD'
colnames(rs3AF)[ncol(rs3AF)]='MeanFD'
colnames(rs4AF)[ncol(rs4AF)]='MeanFD'
colnames(rs5AF)[ncol(rs5AF)]='MeanFD'
colnames(rs6AF)[ncol(rs6AF)]='MeanFD'
colnames(rs1D)[ncol(rs1D)]='MeanFD'
colnames(rs2D)[ncol(rs2D)]='MeanFD'
colnames(rs3D)[ncol(rs3D)]='MeanFD'
colnames(rs4D)[ncol(rs4D)]='MeanFD'
colnames(rs5D)[ncol(rs5D)]='MeanFD'
colnames(rs6D)[ncol(rs6D)]='MeanFD'
# second to last column is RemTRs
colnames(rs1BV)[ncol(rs1BV)-1]='RemTRs'
colnames(rs2BV)[ncol(rs2BV)-1]='RemTRs'
colnames(rs3BV)[ncol(rs3BV)-1]='RemTRs'
colnames(rs4BV)[ncol(rs4BV)-1]='RemTRs'
colnames(rs5BV)[ncol(rs5BV)-1]='RemTRs'
colnames(rs6BV)[ncol(rs6BV)-1]='RemTRs'
colnames(rs1BW)[ncol(rs1BW)-1]='RemTRs'
colnames(rs2BW)[ncol(rs2BW)-1]='RemTRs'
colnames(rs3BW)[ncol(rs3BW)-1]='RemTRs'
colnames(rs4BW)[ncol(rs4BW)-1]='RemTRs'
colnames(rs5BW)[ncol(rs5BW)-1]='RemTRs'
colnames(rs6BW)[ncol(rs6BW)-1]='RemTRs'
colnames(rs1AF)[ncol(rs1AF)-1]='RemTRs'
colnames(rs2AF)[ncol(rs2AF)-1]='RemTRs'
colnames(rs3AF)[ncol(rs3AF)-1]='RemTRs'
colnames(rs4AF)[ncol(rs4AF)-1]='RemTRs'
colnames(rs5AF)[ncol(rs5AF)-1]='RemTRs'
colnames(rs6AF)[ncol(rs6AF)-1]='RemTRs'
colnames(rs1D)[ncol(rs1D)-1]='RemTRs'
colnames(rs2D)[ncol(rs2D)-1]='RemTRs'
colnames(rs3D)[ncol(rs3D)-1]='RemTRs'
colnames(rs4D)[ncol(rs4D)-1]='RemTRs'
colnames(rs5D)[ncol(rs5D)-1]='RemTRs'
colnames(rs6D)[ncol(rs6D)-1]='RemTRs'
# add subject names
subj_order=read.delim('~/rs_Psil_propsMerged_subjOrder_WithSesh.csv',blank.lines.skip = FALSE,sep=',')
# decode drug assignment based on subjSeshDoseCorresp_psilo
# add Task
# code Drug
# combine
# subjects as factor
# Drug as factor
# QC threshold

# Combine all datasets

#################### model each face
for (f in 1:length(
#################### save out stats
