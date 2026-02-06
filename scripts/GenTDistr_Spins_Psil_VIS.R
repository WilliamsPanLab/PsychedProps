# mirror load in from main script, but use it to generate t stats on spun distribution
library(nlme)
numSpins=1900
# initialize t distribution
tdistr_BUP=rep(0,numSpins)
tdistr_Mag=rep(0,numSpins)
tdistr_Seg=rep(0,numSpins)
tdistr_TA=rep(0,numSpins)

# spun props
for (n in 1:numSpins){
	print(n)
	rs1=read.csv(paste0('/scratch/users/apines/data/rs1_Psil_propsMerged_Spin_',n,'_VIS.csv'),header=F)
	rs2=read.csv(paste0('/scratch/users/apines/data/rs2_Psil_propsMerged_Spin_',n,'_VIS.csv'),header=F)
	rs3=read.csv(paste0('/scratch/users/apines/data/rs3_Psil_propsMerged_Spin_',n,'_VIS.csv'),header=F)
	rs4=read.csv(paste0('/scratch/users/apines/data/rs4_Psil_propsMerged_Spin_',n,'_VIS.csv'),header=F)
	rs5=read.csv(paste0('/scratch/users/apines/data/rs5_Psil_propsMerged_Spin_',n,'_VIS.csv'),header=F)
	rs6=read.csv(paste0('/scratch/users/apines/data/rs6_Psil_propsMerged_Spin_',n,'_VIS.csv'),header=F)
	# mimic main analysis code
	rs1$Task='rs'
	rs2$Task='rs2'
	rs3$Task='rs3'
	rs4$Task='rs4'
	rs5$Task='rs5'
	rs6$Task='rs6'
	rs1$OgRow=seq(1:dim(rs1)[1])
	rs2$OgRow=seq(1:dim(rs2)[1])
	rs3$OgRow=seq(1:dim(rs3)[1])
	rs4$OgRow=seq(1:dim(rs4)[1])
	rs5$OgRow=seq(1:dim(rs5)[1])
	rs6$OgRow=seq(1:dim(rs6)[1])
	rs1bv=data.frame(cbind(rs1$V1,rs1$V2,rs1$V3,rs1$V4,rs1$V17,rs1$OgRow,rs1$V21))
	rs1p=data.frame(cbind(rs1$V5,rs1$V6,rs1$V7,rs1$V8,rs1$V18,rs1$OgRow,rs1$V22))
	rs1m1=data.frame(cbind(rs1$V9,rs1$V10,rs1$V11,rs1$V12,rs1$V19,rs1$OgRow,rs1$V23))
	rs1m2=data.frame(cbind(rs1$V13,rs1$V14,rs1$V15,rs1$V16,rs1$V20,rs1$OgRow,rs1$V24))
	colnames(rs1bv)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs','OgRow','FD')
	colnames(rs1p)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs','OgRow','FD')
	colnames(rs1m1)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs','OgRow','FD')
	colnames(rs1m2)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs','OgRow','FD')
	rs2bv=data.frame(cbind(rs2$V1,rs2$V2,rs2$V3,rs2$V4,rs2$V17,rs2$OgRow,rs2$V21))
	rs2p=data.frame(cbind(rs2$V5,rs2$V6,rs2$V7,rs2$V8,rs2$V18,rs2$OgRow,rs2$V22))
	rs2m1=data.frame(cbind(rs2$V9,rs2$V10,rs2$V11,rs2$V12,rs2$V19,rs2$OgRow,rs2$V23))
	rs2m2=data.frame(cbind(rs2$V13,rs2$V14,rs2$V15,rs2$V16,rs2$V20,rs2$OgRow,rs2$V24))
	colnames(rs2bv)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs','OgRow','FD')
	colnames(rs2p)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs','OgRow','FD')
	colnames(rs2m1)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs','OgRow','FD')
	colnames(rs2m2)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs','OgRow','FD')
	rs3bv=data.frame(cbind(rs3$V1,rs3$V2,rs3$V3,rs3$V4,rs3$V17,rs3$OgRow,rs3$V21))
	rs3p=data.frame(cbind(rs3$V5,rs3$V6,rs3$V7,rs3$V8,rs3$V18,rs3$OgRow,rs3$V22))
	rs3m1=data.frame(cbind(rs3$V9,rs3$V10,rs3$V11,rs3$V12,rs3$V19,rs3$OgRow,rs3$V23))
	rs3m2=data.frame(cbind(rs3$V13,rs3$V14,rs3$V15,rs3$V16,rs3$V20,rs3$OgRow,rs3$V24))
	colnames(rs3bv)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs','OgRow','FD')
	colnames(rs3p)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs','OgRow','FD')
	colnames(rs3m1)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs','OgRow','FD')
	colnames(rs3m2)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs','OgRow','FD')
	rs4bv=data.frame(cbind(rs4$V1,rs4$V2,rs4$V3,rs4$V4,rs4$V17,rs4$OgRow,rs4$V21))
	rs4p=data.frame(cbind(rs4$V5,rs4$V6,rs4$V7,rs4$V8,rs4$V18,rs4$OgRow,rs4$V22))
	rs4m1=data.frame(cbind(rs4$V9,rs4$V10,rs4$V11,rs4$V12,rs4$V19,rs4$OgRow,rs4$V23))
	rs4m2=data.frame(cbind(rs4$V13,rs4$V14,rs4$V15,rs4$V16,rs4$V20,rs4$OgRow,rs4$V24))
	colnames(rs4bv)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs','OgRow','FD')
	colnames(rs4p)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs','OgRow','FD')
	colnames(rs4m1)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs','OgRow','FD')
	colnames(rs4m2)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs','OgRow','FD')
	rs5bv=data.frame(cbind(rs5$V1,rs5$V2,rs5$V3,rs5$V4,rs5$V17,rs5$OgRow,rs5$V21))
	rs5p=data.frame(cbind(rs5$V5,rs5$V6,rs5$V7,rs5$V8,rs5$V18,rs5$OgRow,rs5$V22))
	rs5m1=data.frame(cbind(rs5$V9,rs5$V10,rs5$V11,rs5$V12,rs5$V19,rs5$OgRow,rs5$V23))
	rs5m2=data.frame(cbind(rs5$V13,rs5$V14,rs5$V15,rs5$V16,rs5$V20,rs5$OgRow,rs5$V24))
	colnames(rs5bv)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs','OgRow','FD')
	colnames(rs5p)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs','OgRow','FD')
	colnames(rs5m1)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs','OgRow','FD')
	colnames(rs5m2)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs','OgRow','FD')
	rs6bv=data.frame(cbind(rs6$V1,rs6$V2,rs6$V3,rs6$V4,rs6$V17,rs6$OgRow,rs6$V21))
	rs6p=data.frame(cbind(rs6$V5,rs6$V6,rs6$V7,rs6$V8,rs6$V18,rs6$OgRow,rs6$V22))
	rs6m1=data.frame(cbind(rs6$V9,rs6$V10,rs6$V11,rs6$V12,rs6$V19,rs6$OgRow,rs6$V23))
	rs6m2=data.frame(cbind(rs6$V13,rs6$V14,rs6$V15,rs6$V16,rs6$V20,rs6$OgRow,rs6$V24))
	colnames(rs6bv)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs','OgRow','FD')
	colnames(rs6p)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs','OgRow','FD')
	colnames(rs6m1)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs','OgRow','FD')
	colnames(rs6m2)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs','OgRow','FD')
	# get subject IDs subject order csv (should both be identical)
	subjOrder_rs1=read.csv('/oak/stanford/groups/leanew1/users/apines/data/rs1_Psil_propsMerged_subjOrder.csv',blank.lines.skip = FALSE)
	subjOrder_rs2=read.csv('/oak/stanford/groups/leanew1/users/apines/data/rs2_Psil_propsMerged_subjOrder.csv',blank.lines.skip = FALSE)
	subjOrder_rs3=read.csv('/oak/stanford/groups/leanew1/users/apines/data/rs3_Psil_propsMerged_subjOrder.csv',blank.lines.skip = FALSE)
	subjOrder_rs4=read.csv('/oak/stanford/groups/leanew1/users/apines/data/rs4_Psil_propsMerged_subjOrder.csv',blank.lines.skip = FALSE)
	subjOrder_rs5=read.csv('/oak/stanford/groups/leanew1/users/apines/data/rs5_Psil_propsMerged_subjOrder.csv',blank.lines.skip = FALSE)
	subjOrder_rs6=read.csv('/oak/stanford/groups/leanew1/users/apines/data/rs6_Psil_propsMerged_subjOrder.csv',blank.lines.skip = FALSE)
	rs1bv$Subjects=subjOrder_rs1$SubjNameCol_1
	rs2bv$Subjects=subjOrder_rs2$SubjNameCol_1
	rs3bv$Subjects=subjOrder_rs3$SubjNameCol_1
	rs4bv$Subjects=subjOrder_rs4$SubjNameCol_1
	rs5bv$Subjects=subjOrder_rs5$SubjNameCol_1
	rs6bv$Subjects=subjOrder_rs6$SubjNameCol_1
	rs1p$Subjects=subjOrder_rs1$SubjNameCol_1
	rs2p$Subjects=subjOrder_rs2$SubjNameCol_1
	rs3p$Subjects=subjOrder_rs3$SubjNameCol_1
	rs4p$Subjects=subjOrder_rs4$SubjNameCol_1
	rs5p$Subjects=subjOrder_rs5$SubjNameCol_1
	rs6p$Subjects=subjOrder_rs6$SubjNameCol_1
	rs1m1$Subjects=subjOrder_rs1$SubjNameCol_1
	rs2m1$Subjects=subjOrder_rs2$SubjNameCol_1
	rs3m1$Subjects=subjOrder_rs3$SubjNameCol_1
	rs4m1$Subjects=subjOrder_rs4$SubjNameCol_1
	rs5m1$Subjects=subjOrder_rs5$SubjNameCol_1
	rs6m1$Subjects=subjOrder_rs6$SubjNameCol_1
	rs1m2$Subjects=subjOrder_rs1$SubjNameCol_1
	rs2m2$Subjects=subjOrder_rs2$SubjNameCol_1
	rs3m2$Subjects=subjOrder_rs3$SubjNameCol_1
	rs4m2$Subjects=subjOrder_rs4$SubjNameCol_1
	rs5m2$Subjects=subjOrder_rs5$SubjNameCol_1
	rs6m2$Subjects=subjOrder_rs6$SubjNameCol_1
	# add session straight from matlab readout
	subjOrder_rs_ws=read.delim('~/rs_Psil_propsMerged_subjOrder_WithSesh.csv',blank.lines.skip = FALSE,sep=',')
	rs1bv$Session=subjOrder_rs_ws$SubjNameCol_2
	rs2bv$Session=subjOrder_rs_ws$SubjNameCol_2
	rs3bv$Session=subjOrder_rs_ws$SubjNameCol_2
	rs4bv$Session=subjOrder_rs_ws$SubjNameCol_2
	rs5bv$Session=subjOrder_rs_ws$SubjNameCol_2
	rs6bv$Session=subjOrder_rs_ws$SubjNameCol_2
	rs1p$Session=subjOrder_rs_ws$SubjNameCol_2
	rs2p$Session=subjOrder_rs_ws$SubjNameCol_2
	rs3p$Session=subjOrder_rs_ws$SubjNameCol_2
	rs4p$Session=subjOrder_rs_ws$SubjNameCol_2
	rs5p$Session=subjOrder_rs_ws$SubjNameCol_2
	rs6p$Session=subjOrder_rs_ws$SubjNameCol_2
	rs1m1$Session=subjOrder_rs_ws$SubjNameCol_2
	rs2m1$Session=subjOrder_rs_ws$SubjNameCol_2
	rs3m1$Session=subjOrder_rs_ws$SubjNameCol_2
	rs4m1$Session=subjOrder_rs_ws$SubjNameCol_2
	rs5m1$Session=subjOrder_rs_ws$SubjNameCol_2
	rs6m1$Session=subjOrder_rs_ws$SubjNameCol_2
	rs1m2$Session=subjOrder_rs_ws$SubjNameCol_2
	rs2m2$Session=subjOrder_rs_ws$SubjNameCol_2
	rs3m2$Session=subjOrder_rs_ws$SubjNameCol_2
	rs4m2$Session=subjOrder_rs_ws$SubjNameCol_2
	rs5m2$Session=subjOrder_rs_ws$SubjNameCol_2
	rs6m2$Session=subjOrder_rs_ws$SubjNameCol_2
	# add in task (rs to be made equivalent after motion merge)
	rs1bv$Task='rs'
	rs1p$Task='rs'
	rs1m1$Task='rs'
	rs1m2$Task='rs'
	rs2bv$Task='rs2'
	rs2p$Task='rs2'
	rs2m1$Task='rs2'
	rs2m2$Task='rs2'
	rs3bv$Task='rs3'
	rs3p$Task='rs3'
	rs3m1$Task='rs3'
	rs3m2$Task='rs3'
	rs4bv$Task='rs4'
	rs4p$Task='rs4'
	rs4m1$Task='rs4'
	rs4m2$Task='rs4'
	rs5bv$Task='rs5'
	rs5p$Task='rs5'
	rs5m1$Task='rs5'
	rs5m2$Task='rs5'
	rs6bv$Task='rs6'
	rs6p$Task='rs6'
	rs6m1$Task='rs6'
	rs6m2$Task='rs6'
	# add in dosage
	rs1bv$Dosage='none'
	rs1p$Dosage='none'
	rs1m1$Dosage='none'
	rs1m2$Dosage='Drug'
	rs2bv$Dosage='none'
	rs2p$Dosage='none'
	rs2m1$Dosage='none'
	rs2m2$Dosage='Drug'
	rs3bv$Dosage='none'
	rs3p$Dosage='none'
	rs3m1$Dosage='none'
	rs3m2$Dosage='Drug'
	rs4bv$Dosage='none'
	rs4p$Dosage='none'
	rs4m1$Dosage='none'
	rs4m2$Dosage='Drug'
	rs5bv$Dosage='none'
	rs5p$Dosage='none'
	rs5m1$Dosage='none'
	rs5m2$Dosage='Drug'
	rs6bv$Dosage='none'
	rs6p$Dosage='none'
	rs6m1$Dosage='none'
	rs6m2$Dosage='Drug'
	# parse out only existing rows 
	rs1bv=rs1bv[rs1bv$TDProp1>0,]
	rs1p=rs1p[rs1p$TDProp1>0,]
	rs1m1=rs1m1[rs1m1$TDProp1>0,]
	rs1m2=rs1m2[rs1m2$TDProp1>0,]
	rs2bv=rs2bv[rs2bv$TDProp1>0,]
	rs2p=rs2p[rs2p$TDProp1>0,]
	rs2m1=rs2m1[rs2m1$TDProp1>0,]
	rs2m2=rs2m2[rs2m2$TDProp1>0,]
	rs3bv=rs3bv[rs3bv$TDProp1>0,]
	rs3p=rs3p[rs3p$TDProp1>0,]
	rs3m1=rs3m1[rs3m1$TDProp1>0,]
	rs3m2=rs3m2[rs3m2$TDProp1>0,]
	rs4bv=rs4bv[rs4bv$TDProp1>0,]
	rs4p=rs4p[rs4p$TDProp1>0,]
	rs4m1=rs4m1[rs4m1$TDProp1>0,]
	rs4m2=rs4m2[rs4m2$TDProp1>0,]
	rs5bv=rs5bv[rs5bv$TDProp1>0,]
	rs5p=rs5p[rs5p$TDProp1>0,]
	rs5m1=rs5m1[rs5m1$TDProp1>0,]
	rs5m2=rs5m2[rs5m2$TDProp1>0,]
	rs6bv=rs6bv[rs6bv$TDProp1>0,]
	rs6p=rs6p[rs6p$TDProp1>0,]
	rs6m1=rs6m1[rs6m1$TDProp1>0,]
	rs6m2=rs6m2[rs6m2$TDProp1>0,]
	# decode Methylphenidate vs. psilocybin for each PT
	subjDoseCorresp=read.csv('~/subjSeshDoseCorresp_psilo.csv',header=F)
	subjDoseCorresp=data.frame(t(subjDoseCorresp))
	subjDoseCorresp$X2<-as.numeric(subjDoseCorresp$X2)
	# initialize new drug column in all dfs
	rs1bv$Drug <- rep('none', nrow(rs1bv))
	rs1bv$Drug <- factor(rs1bv$Drug, levels = c('none', 'Methyl', 'Psilo','Drug1','Drug2','Drug'))
	rs2bv$Drug=rep('none', nrow(rs2bv))
	rs2bv$Drug<- factor(rs2bv$Drug, levels = c('none', 'Methyl', 'Psilo','Drug1','Drug2','Drug'))
	rs1p$Drug=rep('none', nrow(rs1p))
	rs1p$Drug<- factor(rs1p$Drug, levels = c('none', 'Methyl', 'Psilo','Drug1','Drug2','Drug'))
	rs2p$Drug=rep('none', nrow(rs2p))
	rs2p$Drug<- factor(rs2p$Drug, levels = c('none', 'Methyl', 'Psilo','Drug1','Drug2','Drug'))
	rs1m1$Drug=rep('none', nrow(rs1m1))
	rs1m1$Drug<- factor(rs1m1$Drug, levels = c('none', 'Methyl', 'Psilo','Drug1','Drug2','Drug'))
	rs2m1$Drug=rep('none', nrow(rs2m1))
	rs2m1$Drug<- factor(rs2m1$Drug, levels = c('none', 'Methyl', 'Psilo','Drug1','Drug2','Drug'))
	rs1m2$Drug=rep('none', nrow(rs1m2))
	rs1m2$Drug<- factor(rs1m2$Drug, levels = c('none', 'Methyl', 'Psilo','Drug1','Drug2','Drug'))
	rs2m2$Drug=rep('none', nrow(rs2m2))
	rs2m2$Drug<- factor(rs2m2$Drug, levels = c('none', 'Methyl', 'Psilo','Drug1','Drug2','Drug'))
	rs3bv$Drug=rep('none', nrow(rs3bv))
	rs3bv$Drug<- factor(rs3bv$Drug, levels = c('none', 'Methyl', 'Psilo','Drug1','Drug2','Drug'))
	rs3p$Drug=rep('none', nrow(rs3p))
	rs3p$Drug<- factor(rs3p$Drug, levels = c('none', 'Methyl', 'Psilo','Drug1','Drug2','Drug'))
	rs3m1$Drug=rep('none', nrow(rs3m1))
	rs3m1$Drug<- factor(rs3m1$Drug, levels = c('none', 'Methyl', 'Psilo','Drug1','Drug2','Drug'))
	rs3m2$Drug=rep('none', nrow(rs3m2))
	rs3m2$Drug<- factor(rs3m2$Drug, levels = c('none', 'Methyl', 'Psilo','Drug1','Drug2','Drug'))
	rs4bv$Drug=rep('none', nrow(rs4bv))
	rs4bv$Drug<- factor(rs4bv$Drug, levels = c('none', 'Methyl', 'Psilo','Drug1','Drug2','Drug'))
	rs4p$Drug=rep('none', nrow(rs4p))
	rs4p$Drug<- factor(rs4p$Drug, levels = c('none', 'Methyl', 'Psilo','Drug1','Drug2','Drug'))
	rs4m1$Drug=rep('none', nrow(rs4m1))
	rs4m1$Drug<- factor(rs4m1$Drug, levels = c('none', 'Methyl', 'Psilo','Drug1','Drug2','Drug'))
	rs4m2$Drug=rep('none', nrow(rs4m2))
	rs4m2$Drug<- factor(rs4m2$Drug, levels = c('none', 'Methyl', 'Psilo','Drug1','Drug2','Drug'))
	rs5bv$Drug=rep('none', nrow(rs5bv))
	rs5bv$Drug<- factor(rs5bv$Drug, levels = c('none', 'Methyl', 'Psilo','Drug1','Drug2','Drug'))
	rs5p$Drug=rep('none', nrow(rs5p))
	rs5p$Drug<- factor(rs5p$Drug, levels = c('none', 'Methyl', 'Psilo','Drug1','Drug2','Drug'))
	rs5m1$Drug=rep('none', nrow(rs5m1))
	rs5m1$Drug<- factor(rs5m1$Drug, levels = c('none', 'Methyl', 'Psilo','Drug1','Drug2','Drug'))
	rs5m2$Drug=rep('none', nrow(rs5m2))
	rs5m2$Drug<- factor(rs5m2$Drug, levels = c('none', 'Methyl', 'Psilo','Drug1','Drug2','Drug'))
	rs6bv$Drug=rep('none', nrow(rs6bv))
	rs6bv$Drug<- factor(rs6bv$Drug, levels = c('none', 'Methyl', 'Psilo','Drug1','Drug2','Drug'))
	rs6p$Drug=rep('none', nrow(rs6p))
	rs6p$Drug<- factor(rs6p$Drug, levels = c('none', 'Methyl', 'Psilo','Drug1','Drug2','Drug'))
	rs6m1$Drug=rep('none', nrow(rs6m1))
	rs6m1$Drug<- factor(rs6m1$Drug, levels = c('none', 'Methyl', 'Psilo','Drug1','Drug2','Drug'))
	rs6m2$Drug=rep('none', nrow(rs6m2))
	rs6m2$Drug<- factor(rs6m2$Drug, levels = c('none', 'Methyl', 'Psilo','Drug1','Drug2','Drug'))
	# for resting-state 1
	for (s in 1:length(unique(rs1m2$Subjects))){
		  # this subject
		  subj=unique(rs1m2$Subjects)[s]
	    # Subset rows for the current subject
	    subjRows <- rs1m2[rs1m2$Subjects == subj, ]
	    # Subset rows for the 'Drug' dosage
	    subjDrugRows <- subjRows[subjRows$Dosage == 'Drug', ]
	      # Identify 'Drug1' and 'Drug2' based on 'OgRow'
	      drug1Rows <- which(subjDrugRows$Session == 1)
	      drug2Rows <- which(subjDrugRows$Session == 2)
	          subjRows$Drug[drug1Rows] <- 'Drug1'
	          subjRows$Drug[drug2Rows] <- 'Drug2'
		    ### Now pull out methyl vs. psilo
		    DoseCorresp=subjDoseCorresp[subjDoseCorresp$X1==subj,]
		    # if Drug1 = 1 in key, Drug1 = psilo
		    if (subjDrugRows$Drug[1] == 'DrugOnly_CantSeeMe') {
			        subjRows$Drug[subjRows$Drug == 'Drug'] <- 'Psilo'
		      } else if (DoseCorresp$X2 == 1) {
			          subjRows$Drug[subjRows$Drug == 'Drug1'] <- 'Psilo'
		          subjRows$Drug[subjRows$Drug == 'Drug2'] <- 'Methyl'
			    } else if (DoseCorresp$X2 == 2) {
				        subjRows$Drug[subjRows$Drug == 'Drug1'] <- 'Methyl'
			      subjRows$Drug[subjRows$Drug == 'Drug2'] <- 'Psilo'
			        } 
		      # find where OG row corresponds to methyl
		      MethylRows=subjRows$OgRow[subjRows$Drug=='Methyl']
		      # and psilo
		      PsiloRows=subjRows$OgRow[subjRows$Drug=='Psilo']
		        # set
		        rs1m2$Drug[rs1m2$OgRow==MethylRows[1]]='Methyl'
		        rs1m2$Drug[rs1m2$OgRow==PsiloRows[1]]='Psilo'
	}
	# for resting-state 2
	for (s in 1:length(unique(rs2m2$Subjects))){
		  # this subject
		  subj=unique(rs2m2$Subjects)[s]
	    # Subset rows for the current subject
	    subjRows <- rs2m2[rs2m2$Subjects == subj, ]
	    # Subset rows for the 'Drug' dosage
	    subjDrugRows <- subjRows[subjRows$Dosage == 'Drug', ]
	      # Identify 'Drug1' and 'Drug2' based on 'OgRow'
	      drug1Rows <- which(subjDrugRows$Session == 1)
	      drug2Rows <- which(subjDrugRows$Session == 2)
	        #  # Assign 'Drug' values to the original rows
	          subjRows$Drug[drug1Rows] <- 'Drug1'
	          subjRows$Drug[drug2Rows] <- 'Drug2'
		    ### Now pull out methyl vs. psilo
		    DoseCorresp=subjDoseCorresp[subjDoseCorresp$X1==subj,]
		    # if Drug1 = 1 in key, Drug1 = psilo
		    if (subjDrugRows$Drug[1] == 'DrugOnly_CantSeeMe') {
			        subjRows$Drug[subjRows$Drug == 'Drug'] <- 'Psilo'
		      } else if (DoseCorresp$X2 == 1) {
			          subjRows$Drug[subjRows$Drug == 'Drug1'] <- 'Psilo'
		          subjRows$Drug[subjRows$Drug == 'Drug2'] <- 'Methyl'
			    } else if (DoseCorresp$X2 == 2) {
				        subjRows$Drug[subjRows$Drug == 'Drug1'] <- 'Methyl'
			      subjRows$Drug[subjRows$Drug == 'Drug2'] <- 'Psilo'
			        } 
		      # find where OG row corresponds to methyl
		      MethylRows=subjRows$OgRow[subjRows$Drug=='Methyl']
		      # and psilo
		      PsiloRows=subjRows$OgRow[subjRows$Drug=='Psilo']
		        # set
		        rs2m2$Drug[rs2m2$OgRow==MethylRows[1]]='Methyl'
		        rs2m2$Drug[rs2m2$OgRow==PsiloRows[1]]='Psilo'
	}
	# for resting-state 3
	for (s in 1:length(unique(rs3m2$Subjects))){
		  # this subject
		  subj=unique(rs3m2$Subjects)[s]
	    # Subset rows for the current subject
	    subjRows <- rs3m2[rs3m2$Subjects == subj, ]
	    # Subset rows for the 'Drug' dosage
	    subjDrugRows <- subjRows[subjRows$Dosage == 'Drug', ]
	        # Identify 'Drug1' and 'Drug2' based on 'OgRow'
	      drug1Rows <- which(subjDrugRows$Session == 1)
	      drug2Rows <- which(subjDrugRows$Session == 2)
	        # if min=max, label only 1 drug
	        #  # Assign 'Drug' values to the original rows
	          subjRows$Drug[drug1Rows] <- 'Drug1'
	          subjRows$Drug[drug2Rows] <- 'Drug2'
		    #}
		    ### Now pull out methyl vs. psilo
		    DoseCorresp=subjDoseCorresp[subjDoseCorresp$X1==subj,]
		    # if Drug1 = 1 in key, Drug1 = psilo
		    if (subjDrugRows$Drug[1] == 'DrugOnly_cantSeeMe') {
			        subjRows$Drug[subjRows$Drug == 'Drug'] <- 'Psilo'
		      } else if (DoseCorresp$X2 == 1) {
			          subjRows$Drug[subjRows$Drug == 'Drug1'] <- 'Psilo'
		          subjRows$Drug[subjRows$Drug == 'Drug2'] <- 'Methyl'
			    } else if (DoseCorresp$X2 == 2) {
				        subjRows$Drug[subjRows$Drug == 'Drug1'] <- 'Methyl'
			      subjRows$Drug[subjRows$Drug == 'Drug2'] <- 'Psilo'
			        }
		      # find where OG row corresponds to methyl
		      MethylRows=subjRows$OgRow[subjRows$Drug=='Methyl']
		      # and psilo
		      PsiloRows=subjRows$OgRow[subjRows$Drug=='Psilo']
		        # set
		        rs3m2$Drug[rs3m2$OgRow==MethylRows[1]]='Methyl'
		        rs3m2$Drug[rs3m2$OgRow==PsiloRows[1]]='Psilo'
	}
	# for resting-state 4
	for (s in 1:length(unique(rs4m2$Subjects))){
		  # this subject
		  subj=unique(rs4m2$Subjects)[s]
	    # Subset rows for the current subject
	    subjRows <- rs4m2[rs4m2$Subjects == subj, ]
	    # Subset rows for the 'Drug' dosage
	    subjDrugRows <- subjRows[subjRows$Dosage == 'Drug', ]
	      # Identify 'Drug1' and 'Drug2' based on 'OgRow'
	      drug1Rows <- which(subjDrugRows$Session == 1)
	      drug2Rows <- which(subjDrugRows$Session == 2)
	        # if min=max, label only 1 drug
	        #  # Assign 'Drug' values to the original rows
	          subjRows$Drug[drug1Rows] <- 'Drug1'
	          subjRows$Drug[drug2Rows] <- 'Drug2'
		    #}
		    ### Now pull out methyl vs. psilo
		    DoseCorresp=subjDoseCorresp[subjDoseCorresp$X1==subj,]
		    # if Drug1 = 1 in key, Drug1 = psilo
		    if (subjDrugRows$Drug[1] == 'DrugOnly_CantSeeMe') {
			        subjRows$Drug[subjRows$Drug == 'Drug'] <- 'Psilo'
		      } else if (DoseCorresp$X2 == 1) {
			          subjRows$Drug[subjRows$Drug == 'Drug1'] <- 'Psilo'
		          subjRows$Drug[subjRows$Drug == 'Drug2'] <- 'Methyl'
			    } else if (DoseCorresp$X2 == 2) {
				        subjRows$Drug[subjRows$Drug == 'Drug1'] <- 'Methyl'
			      subjRows$Drug[subjRows$Drug == 'Drug2'] <- 'Psilo'
			        }
		      # find where OG row corresponds to methyl
		      MethylRows=subjRows$OgRow[subjRows$Drug=='Methyl']
		      # and psilo
		      PsiloRows=subjRows$OgRow[subjRows$Drug=='Psilo']
		        # set
		        rs4m2$Drug[rs4m2$OgRow==MethylRows[1]]='Methyl'
		        rs4m2$Drug[rs4m2$OgRow==PsiloRows[1]]='Psilo'
	}
	# for resting-state 5
	for (s in 1:length(unique(rs5m2$Subjects))){
		  # this subject
		  subj=unique(rs5m2$Subjects)[s]
	    # Subset rows for the current subject
	    subjRows <- rs5m2[rs5m2$Subjects == subj, ]
	    # Subset rows for the 'Drug' dosage
	    subjDrugRows <- subjRows[subjRows$Dosage == 'Drug', ]
	      # Identify 'Drug1' and 'Drug2' based on 'OgRow'
	      drug1Rows <- which(subjDrugRows$Session == 1)
	      drug2Rows <- which(subjDrugRows$Session == 2)
	        #  # Assign 'Drug' values to the original rows
	          subjRows$Drug[drug1Rows] <- 'Drug1'
	          subjRows$Drug[drug2Rows] <- 'Drug2'
		    #}
		    ### Now pull out methyl vs. psilo
		    DoseCorresp=subjDoseCorresp[subjDoseCorresp$X1==subj,]
		    # if Drug1 = 1 in key, Drug1 = psilo
		    if (subjDrugRows$Drug[1] == 'DrugOnly_CantSeeMe') {
			        subjRows$Drug[subjRows$Drug == 'Drug'] <- 'Psilo'
		      } else if (DoseCorresp$X2 == 1) {
			          subjRows$Drug[subjRows$Drug == 'Drug1'] <- 'Psilo'
		          subjRows$Drug[subjRows$Drug == 'Drug2'] <- 'Methyl'
			    } else if (DoseCorresp$X2 == 2) {
				        subjRows$Drug[subjRows$Drug == 'Drug1'] <- 'Methyl'
			      subjRows$Drug[subjRows$Drug == 'Drug2'] <- 'Psilo'
			        }
		      # find where OG row corresponds to methyl
		      MethylRows=subjRows$OgRow[subjRows$Drug=='Methyl']
		      # and psilo
		      PsiloRows=subjRows$OgRow[subjRows$Drug=='Psilo']
		        # set
		        rs5m2$Drug[rs5m2$OgRow==MethylRows[1]]='Methyl'
		        rs5m2$Drug[rs5m2$OgRow==PsiloRows[1]]='Psilo'
	}
	# for resting-state 6
	for (s in 1:length(unique(rs6m2$Subjects))){
		  # this subject
		  subj=unique(rs6m2$Subjects)[s]
	    # Subset rows for the current subject
	    subjRows <- rs6m2[rs6m2$Subjects == subj, ]
	    # Subset rows for the 'Drug' dosage
	    subjDrugRows <- subjRows[subjRows$Dosage == 'Drug', ]
	      # Identify 'Drug1' and 'Drug2' based on 'OgRow'
	      drug1Rows <- which(subjDrugRows$Session == 1)
	      drug2Rows <- which(subjDrugRows$Session == 2)
	        #  # Assign 'Drug' values to the original rows
	          subjRows$Drug[drug1Rows] <- 'Drug1'
	          subjRows$Drug[drug2Rows] <- 'Drug2'
		    ### Now pull out methyl vs. psilo
		    DoseCorresp=subjDoseCorresp[subjDoseCorresp$X1==subj,]
		    # if Drug1 = 1 in key, Drug1 = psilo
		    if (subjDrugRows$Drug[1] == 'DrugOnly_CantSeeMe') {
			        subjRows$Drug[subjRows$Drug == 'Drug'] <- 'Psilo'
		      } else if (DoseCorresp$X2 == 1) {
			          subjRows$Drug[subjRows$Drug == 'Drug1'] <- 'Psilo'
		          subjRows$Drug[subjRows$Drug == 'Drug2'] <- 'Methyl'
			    } else if (DoseCorresp$X2 == 2) {
				        subjRows$Drug[subjRows$Drug == 'Drug1'] <- 'Methyl'
			      subjRows$Drug[subjRows$Drug == 'Drug2'] <- 'Psilo'
			        }
		      # find where OG row corresponds to methyl
		      MethylRows=subjRows$OgRow[subjRows$Drug=='Methyl']
		      # and psilo
		      PsiloRows=subjRows$OgRow[subjRows$Drug=='Psilo']
		        # set
		        rs6m2$Drug[rs6m2$OgRow==MethylRows[1]]='Methyl'
		        rs6m2$Drug[rs6m2$OgRow==PsiloRows[1]]='Psilo'
	}

	# combine all
	allScans=rbind(rs1bv,rs1p,rs1m1,rs1m2,rs2bv,rs2p,rs2m1,rs2m2,rs3bv,rs3p,rs3m1,rs3m2,rs4bv,rs4p,rs4m1,rs4m2,rs5bv,rs5p,rs5m1,rs5m2,rs6bv,rs6p,rs6m1,rs6m2)
	# before scans
	Before=rbind(rs1bv,rs2bv,rs3bv,rs4bv,rs5bv,rs6bv)
	Before$Chronology='Before'
	# after
	After=rbind(rs1m1,rs2m1,rs3m1,rs4m1,rs5m1,rs6m1)
	After$Chronology='After'
	B_A=rbind(Before,After)
	# set non-drug conditions to none
	allScans$Drug[allScans$Dosage=='none']='none'
	# retain a legacy subjects naming convention for MEQ use later
	allScans$LegacySubjs=allScans$Subjects
	# final thing is to convert all "before" to "Baseline", "between" to "Between". and "after" to "After"
	allScans$Session <- gsub("before", "Baseline", allScans$Session)
	allScans$Session <- gsub("between", "Between", allScans$Session)
	allScans$Session <- gsub("after", "After", allScans$Session)
	# qc criteria
	allScans=allScans[allScans$RemTRs>250,]
	# model
	allScansNoMeth=allScans[allScans$Drug!='Methyl',]
	fit_lme <- lme(TDProp1 ~ Drug + RemTRs + FD, random = ~ 1 | Subjects, data = allScansNoMeth)
	summaryLME<-summary(fit_lme)
	tdistr_BUP[n]=summaryLME$tTable['DrugPsilo','t-value']
}
print('done with spun angles')
######################
# magnitudes
################
for (n in 1:numSpins){
 	rs1=read.csv(paste0('/scratch/users/apines/data/rs1_Psil_VISMagMerged_Spin_',n,'.csv'),header=F)
        rs2=read.csv(paste0('/scratch/users/apines/data/rs2_Psil_VISMagMerged_Spin_',n,'.csv'),header=F)
        rs3=read.csv(paste0('/scratch/users/apines/data/rs3_Psil_VISMagMerged_Spin_',n,'.csv'),header=F)
        rs4=read.csv(paste0('/scratch/users/apines/data/rs4_Psil_VISMagMerged_Spin_',n,'.csv'),header=F)
        rs5=read.csv(paste0('/scratch/users/apines/data/rs5_Psil_VISMagMerged_Spin_',n,'.csv'),header=F)
        rs6=read.csv(paste0('/scratch/users/apines/data/rs6_Psil_VISMagMerged_Spin_',n,'.csv'),header=F)
        # mimic main analysis code
        rs1$Task='rs'
        rs2$Task='rs2'
        rs3$Task='rs3'
        rs4$Task='rs4'
        rs5$Task='rs5'
        rs6$Task='rs6'
        rs1$OgRow=seq(1:dim(rs1)[1])
        rs2$OgRow=seq(1:dim(rs2)[1])
        rs3$OgRow=seq(1:dim(rs3)[1])
        rs4$OgRow=seq(1:dim(rs4)[1])
        rs5$OgRow=seq(1:dim(rs5)[1])
        rs6$OgRow=seq(1:dim(rs6)[1])
        rs1bv=data.frame(cbind(rs1$V1,rs1$V2,rs1$V3,rs1$V4,rs1$V17,rs1$OgRow,rs1$V21))
        rs1p=data.frame(cbind(rs1$V5,rs1$V6,rs1$V7,rs1$V8,rs1$V18,rs1$OgRow,rs1$V22))
        rs1m1=data.frame(cbind(rs1$V9,rs1$V10,rs1$V11,rs1$V12,rs1$V19,rs1$OgRow,rs1$V23))
        rs1m2=data.frame(cbind(rs1$V13,rs1$V14,rs1$V15,rs1$V16,rs1$V20,rs1$OgRow,rs1$V24))
        colnames(rs1bv)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs','OgRow','FD')
        colnames(rs1p)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs','OgRow','FD')
        colnames(rs1m1)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs','OgRow','FD')
        colnames(rs1m2)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs','OgRow','FD')
        rs2bv=data.frame(cbind(rs2$V1,rs2$V2,rs2$V3,rs2$V4,rs2$V17,rs2$OgRow,rs2$V21))
        rs2p=data.frame(cbind(rs2$V5,rs2$V6,rs2$V7,rs2$V8,rs2$V18,rs2$OgRow,rs2$V22))
        rs2m1=data.frame(cbind(rs2$V9,rs2$V10,rs2$V11,rs2$V12,rs2$V19,rs2$OgRow,rs2$V23))
        rs2m2=data.frame(cbind(rs2$V13,rs2$V14,rs2$V15,rs2$V16,rs2$V20,rs2$OgRow,rs2$V24))
        colnames(rs2bv)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs','OgRow','FD')
        colnames(rs2p)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs','OgRow','FD')
        colnames(rs2m1)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs','OgRow','FD')
        colnames(rs2m2)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs','OgRow','FD')
        rs3bv=data.frame(cbind(rs3$V1,rs3$V2,rs3$V3,rs3$V4,rs3$V17,rs3$OgRow,rs3$V21))
        rs3p=data.frame(cbind(rs3$V5,rs3$V6,rs3$V7,rs3$V8,rs3$V18,rs3$OgRow,rs3$V22))
        rs3m1=data.frame(cbind(rs3$V9,rs3$V10,rs3$V11,rs3$V12,rs3$V19,rs3$OgRow,rs3$V23))
        rs3m2=data.frame(cbind(rs3$V13,rs3$V14,rs3$V15,rs3$V16,rs3$V20,rs3$OgRow,rs3$V24))
        colnames(rs3bv)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs','OgRow','FD')
        colnames(rs3p)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs','OgRow','FD')
        colnames(rs3m1)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs','OgRow','FD')
        colnames(rs3m2)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs','OgRow','FD')
        rs4bv=data.frame(cbind(rs4$V1,rs4$V2,rs4$V3,rs4$V4,rs4$V17,rs4$OgRow,rs4$V21))
        rs4p=data.frame(cbind(rs4$V5,rs4$V6,rs4$V7,rs4$V8,rs4$V18,rs4$OgRow,rs4$V22))
        rs4m1=data.frame(cbind(rs4$V9,rs4$V10,rs4$V11,rs4$V12,rs4$V19,rs4$OgRow,rs4$V23))
        rs4m2=data.frame(cbind(rs4$V13,rs4$V14,rs4$V15,rs4$V16,rs4$V20,rs4$OgRow,rs4$V24))
        colnames(rs4bv)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs','OgRow','FD')
        colnames(rs4p)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs','OgRow','FD')
        colnames(rs4m1)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs','OgRow','FD')
        colnames(rs4m2)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs','OgRow','FD')
        rs5bv=data.frame(cbind(rs5$V1,rs5$V2,rs5$V3,rs5$V4,rs5$V17,rs5$OgRow,rs5$V21))
        rs5p=data.frame(cbind(rs5$V5,rs5$V6,rs5$V7,rs5$V8,rs5$V18,rs5$OgRow,rs5$V22))
        rs5m1=data.frame(cbind(rs5$V9,rs5$V10,rs5$V11,rs5$V12,rs5$V19,rs5$OgRow,rs5$V23))
        rs5m2=data.frame(cbind(rs5$V13,rs5$V14,rs5$V15,rs5$V16,rs5$V20,rs5$OgRow,rs5$V24))
        colnames(rs5bv)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs','OgRow','FD')
        colnames(rs5p)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs','OgRow','FD')
        colnames(rs5m1)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs','OgRow','FD')
        colnames(rs5m2)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs','OgRow','FD')
        rs6bv=data.frame(cbind(rs6$V1,rs6$V2,rs6$V3,rs6$V4,rs6$V17,rs6$OgRow,rs6$V21))
        rs6p=data.frame(cbind(rs6$V5,rs6$V6,rs6$V7,rs6$V8,rs6$V18,rs6$OgRow,rs6$V22))
        rs6m1=data.frame(cbind(rs6$V9,rs6$V10,rs6$V11,rs6$V12,rs6$V19,rs6$OgRow,rs6$V23))
        rs6m2=data.frame(cbind(rs6$V13,rs6$V14,rs6$V15,rs6$V16,rs6$V20,rs6$OgRow,rs6$V24))
        colnames(rs6bv)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs','OgRow','FD')
        colnames(rs6p)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs','OgRow','FD')
        colnames(rs6m1)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs','OgRow','FD')
        colnames(rs6m2)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs','OgRow','FD')
        # get subject IDs subject order csv (should both be identical)
        subjOrder_rs1=read.csv('/oak/stanford/groups/leanew1/users/apines/data/rs1_Psil_propsMerged_subjOrder.csv',blank.lines.skip = FALSE)
        subjOrder_rs2=read.csv('/oak/stanford/groups/leanew1/users/apines/data/rs2_Psil_propsMerged_subjOrder.csv',blank.lines.skip = FALSE)
        subjOrder_rs3=read.csv('/oak/stanford/groups/leanew1/users/apines/data/rs3_Psil_propsMerged_subjOrder.csv',blank.lines.skip = FALSE)
        subjOrder_rs4=read.csv('/oak/stanford/groups/leanew1/users/apines/data/rs4_Psil_propsMerged_subjOrder.csv',blank.lines.skip = FALSE)
        subjOrder_rs5=read.csv('/oak/stanford/groups/leanew1/users/apines/data/rs5_Psil_propsMerged_subjOrder.csv',blank.lines.skip = FALSE)
        subjOrder_rs6=read.csv('/oak/stanford/groups/leanew1/users/apines/data/rs6_Psil_propsMerged_subjOrder.csv',blank.lines.skip = FALSE)
        rs1bv$Subjects=subjOrder_rs1$SubjNameCol_1
        rs2bv$Subjects=subjOrder_rs2$SubjNameCol_1
        rs3bv$Subjects=subjOrder_rs3$SubjNameCol_1
        rs4bv$Subjects=subjOrder_rs4$SubjNameCol_1
        rs5bv$Subjects=subjOrder_rs5$SubjNameCol_1
        rs6bv$Subjects=subjOrder_rs6$SubjNameCol_1
        rs1p$Subjects=subjOrder_rs1$SubjNameCol_1
        rs2p$Subjects=subjOrder_rs2$SubjNameCol_1
        rs3p$Subjects=subjOrder_rs3$SubjNameCol_1
        rs4p$Subjects=subjOrder_rs4$SubjNameCol_1
        rs5p$Subjects=subjOrder_rs5$SubjNameCol_1
        rs6p$Subjects=subjOrder_rs6$SubjNameCol_1
        rs1m1$Subjects=subjOrder_rs1$SubjNameCol_1
        rs2m1$Subjects=subjOrder_rs2$SubjNameCol_1
        rs3m1$Subjects=subjOrder_rs3$SubjNameCol_1
        rs4m1$Subjects=subjOrder_rs4$SubjNameCol_1
        rs5m1$Subjects=subjOrder_rs5$SubjNameCol_1
        rs6m1$Subjects=subjOrder_rs6$SubjNameCol_1
        rs1m2$Subjects=subjOrder_rs1$SubjNameCol_1
        rs2m2$Subjects=subjOrder_rs2$SubjNameCol_1
        rs3m2$Subjects=subjOrder_rs3$SubjNameCol_1
        rs4m2$Subjects=subjOrder_rs4$SubjNameCol_1
        rs5m2$Subjects=subjOrder_rs5$SubjNameCol_1
        rs6m2$Subjects=subjOrder_rs6$SubjNameCol_1
        # add session straight from matlab readout
        subjOrder_rs_ws=read.delim('~/rs_Psil_propsMerged_subjOrder_WithSesh.csv',blank.lines.skip = FALSE,sep=',')
        rs1bv$Session=subjOrder_rs_ws$SubjNameCol_2
        rs2bv$Session=subjOrder_rs_ws$SubjNameCol_2
        rs3bv$Session=subjOrder_rs_ws$SubjNameCol_2
        rs4bv$Session=subjOrder_rs_ws$SubjNameCol_2
        rs5bv$Session=subjOrder_rs_ws$SubjNameCol_2
        rs6bv$Session=subjOrder_rs_ws$SubjNameCol_2
        rs1p$Session=subjOrder_rs_ws$SubjNameCol_2
        rs2p$Session=subjOrder_rs_ws$SubjNameCol_2
        rs3p$Session=subjOrder_rs_ws$SubjNameCol_2
        rs4p$Session=subjOrder_rs_ws$SubjNameCol_2
        rs5p$Session=subjOrder_rs_ws$SubjNameCol_2
        rs6p$Session=subjOrder_rs_ws$SubjNameCol_2
        rs1m1$Session=subjOrder_rs_ws$SubjNameCol_2
        rs2m1$Session=subjOrder_rs_ws$SubjNameCol_2
        rs3m1$Session=subjOrder_rs_ws$SubjNameCol_2
        rs4m1$Session=subjOrder_rs_ws$SubjNameCol_2
        rs5m1$Session=subjOrder_rs_ws$SubjNameCol_2
        rs6m1$Session=subjOrder_rs_ws$SubjNameCol_2
        rs1m2$Session=subjOrder_rs_ws$SubjNameCol_2
        rs2m2$Session=subjOrder_rs_ws$SubjNameCol_2
        rs3m2$Session=subjOrder_rs_ws$SubjNameCol_2
        rs4m2$Session=subjOrder_rs_ws$SubjNameCol_2
        rs5m2$Session=subjOrder_rs_ws$SubjNameCol_2
        rs6m2$Session=subjOrder_rs_ws$SubjNameCol_2
        # add in task (rs to be made equivalent after motion merge)
        rs1bv$Task='rs'
        rs1p$Task='rs'
        rs1m1$Task='rs'
        rs1m2$Task='rs'
        rs2bv$Task='rs2'
        rs2p$Task='rs2'
        rs2m1$Task='rs2'
        rs2m2$Task='rs2'
        rs3bv$Task='rs3'
        rs3p$Task='rs3'
        rs3m1$Task='rs3'
        rs3m2$Task='rs3'
        rs4bv$Task='rs4'
        rs4p$Task='rs4'
        rs4m1$Task='rs4'
        rs4m2$Task='rs4'
        rs5bv$Task='rs5'
        rs5p$Task='rs5'
        rs5m1$Task='rs5'
        rs5m2$Task='rs5'
        rs6bv$Task='rs6'
        rs6p$Task='rs6'
        rs6m1$Task='rs6'
        rs6m2$Task='rs6'
        # add in dosage
        rs1bv$Dosage='none'
        rs1p$Dosage='none'
        rs1m1$Dosage='none'
        rs1m2$Dosage='Drug'
        rs2bv$Dosage='none'
        rs2p$Dosage='none'
        rs2m1$Dosage='none'
        rs2m2$Dosage='Drug'
        rs3bv$Dosage='none'
        rs3p$Dosage='none'
        rs3m1$Dosage='none'
        rs3m2$Dosage='Drug'
        rs4bv$Dosage='none'
        rs4p$Dosage='none'
        rs4m1$Dosage='none'
        rs4m2$Dosage='Drug'
        rs5bv$Dosage='none'
        rs5p$Dosage='none'
        rs5m1$Dosage='none'
        rs5m2$Dosage='Drug'
        rs6bv$Dosage='none'
        rs6p$Dosage='none'
        rs6m1$Dosage='none'
        rs6m2$Dosage='Drug'
        # parse out only existing rows
        rs1bv=rs1bv[rs1bv$TDProp1>0,]
        rs1p=rs1p[rs1p$TDProp1>0,]
        rs1m1=rs1m1[rs1m1$TDProp1>0,]
        rs1m2=rs1m2[rs1m2$TDProp1>0,]
        rs2bv=rs2bv[rs2bv$TDProp1>0,]
        rs2p=rs2p[rs2p$TDProp1>0,]
        rs2m1=rs2m1[rs2m1$TDProp1>0,]
        rs2m2=rs2m2[rs2m2$TDProp1>0,]
        rs3bv=rs3bv[rs3bv$TDProp1>0,]
        rs3p=rs3p[rs3p$TDProp1>0,]
        rs3m1=rs3m1[rs3m1$TDProp1>0,]
        rs3m2=rs3m2[rs3m2$TDProp1>0,]
        rs4bv=rs4bv[rs4bv$TDProp1>0,]
        rs4p=rs4p[rs4p$TDProp1>0,]
        rs4m1=rs4m1[rs4m1$TDProp1>0,]
        rs4m2=rs4m2[rs4m2$TDProp1>0,]
        rs5bv=rs5bv[rs5bv$TDProp1>0,]
        rs5p=rs5p[rs5p$TDProp1>0,]
        rs5m1=rs5m1[rs5m1$TDProp1>0,]
        rs5m2=rs5m2[rs5m2$TDProp1>0,]
        rs6bv=rs6bv[rs6bv$TDProp1>0,]
        rs6p=rs6p[rs6p$TDProp1>0,]
        rs6m1=rs6m1[rs6m1$TDProp1>0,]
        rs6m2=rs6m2[rs6m2$TDProp1>0,]
        # decode Methylphenidate vs. psilocybin for each PT
        subjDoseCorresp=read.csv('~/subjSeshDoseCorresp_psilo.csv',header=F)
        subjDoseCorresp=data.frame(t(subjDoseCorresp))
        subjDoseCorresp$X2<-as.numeric(subjDoseCorresp$X2)
        # initialize new drug column in all dfs
        rs1bv$Drug <- rep('none', nrow(rs1bv))
        rs1bv$Drug <- factor(rs1bv$Drug, levels = c('none', 'Methyl', 'Psilo','Drug1','Drug2','Drug'))
        rs2bv$Drug=rep('none', nrow(rs2bv))
        rs2bv$Drug<- factor(rs2bv$Drug, levels = c('none', 'Methyl', 'Psilo','Drug1','Drug2','Drug'))
        rs1p$Drug=rep('none', nrow(rs1p))
        rs1p$Drug<- factor(rs1p$Drug, levels = c('none', 'Methyl', 'Psilo','Drug1','Drug2','Drug'))
        rs2p$Drug=rep('none', nrow(rs2p))
        rs2p$Drug<- factor(rs2p$Drug, levels = c('none', 'Methyl', 'Psilo','Drug1','Drug2','Drug'))
        rs1m1$Drug=rep('none', nrow(rs1m1))
        rs1m1$Drug<- factor(rs1m1$Drug, levels = c('none', 'Methyl', 'Psilo','Drug1','Drug2','Drug'))
        rs2m1$Drug=rep('none', nrow(rs2m1))
        rs2m1$Drug<- factor(rs2m1$Drug, levels = c('none', 'Methyl', 'Psilo','Drug1','Drug2','Drug'))
        rs1m2$Drug=rep('none', nrow(rs1m2))
        rs1m2$Drug<- factor(rs1m2$Drug, levels = c('none', 'Methyl', 'Psilo','Drug1','Drug2','Drug'))
        rs2m2$Drug=rep('none', nrow(rs2m2))
        rs2m2$Drug<- factor(rs2m2$Drug, levels = c('none', 'Methyl', 'Psilo','Drug1','Drug2','Drug'))
        rs3bv$Drug=rep('none', nrow(rs3bv))
        rs3bv$Drug<- factor(rs3bv$Drug, levels = c('none', 'Methyl', 'Psilo','Drug1','Drug2','Drug'))
        rs3p$Drug=rep('none', nrow(rs3p))
        rs3p$Drug<- factor(rs3p$Drug, levels = c('none', 'Methyl', 'Psilo','Drug1','Drug2','Drug'))
        rs3m1$Drug=rep('none', nrow(rs3m1))
        rs3m1$Drug<- factor(rs3m1$Drug, levels = c('none', 'Methyl', 'Psilo','Drug1','Drug2','Drug'))
        rs3m2$Drug=rep('none', nrow(rs3m2))
        rs3m2$Drug<- factor(rs3m2$Drug, levels = c('none', 'Methyl', 'Psilo','Drug1','Drug2','Drug'))
        rs4bv$Drug=rep('none', nrow(rs4bv))
        rs4bv$Drug<- factor(rs4bv$Drug, levels = c('none', 'Methyl', 'Psilo','Drug1','Drug2','Drug'))
        rs4p$Drug=rep('none', nrow(rs4p))
        rs4p$Drug<- factor(rs4p$Drug, levels = c('none', 'Methyl', 'Psilo','Drug1','Drug2','Drug'))
        rs4m1$Drug=rep('none', nrow(rs4m1))
        rs4m1$Drug<- factor(rs4m1$Drug, levels = c('none', 'Methyl', 'Psilo','Drug1','Drug2','Drug'))
        rs4m2$Drug=rep('none', nrow(rs4m2))
        rs4m2$Drug<- factor(rs4m2$Drug, levels = c('none', 'Methyl', 'Psilo','Drug1','Drug2','Drug'))
        rs5bv$Drug=rep('none', nrow(rs5bv))
        rs5bv$Drug<- factor(rs5bv$Drug, levels = c('none', 'Methyl', 'Psilo','Drug1','Drug2','Drug'))
        rs5p$Drug=rep('none', nrow(rs5p))
        rs5p$Drug<- factor(rs5p$Drug, levels = c('none', 'Methyl', 'Psilo','Drug1','Drug2','Drug'))
        rs5m1$Drug=rep('none', nrow(rs5m1))
        rs5m1$Drug<- factor(rs5m1$Drug, levels = c('none', 'Methyl', 'Psilo','Drug1','Drug2','Drug'))
        rs5m2$Drug=rep('none', nrow(rs5m2))
        rs5m2$Drug<- factor(rs5m2$Drug, levels = c('none', 'Methyl', 'Psilo','Drug1','Drug2','Drug'))
        rs6bv$Drug=rep('none', nrow(rs6bv))
        rs6bv$Drug<- factor(rs6bv$Drug, levels = c('none', 'Methyl', 'Psilo','Drug1','Drug2','Drug'))
        rs6p$Drug=rep('none', nrow(rs6p))
        rs6p$Drug<- factor(rs6p$Drug, levels = c('none', 'Methyl', 'Psilo','Drug1','Drug2','Drug'))
        rs6m1$Drug=rep('none', nrow(rs6m1))
        rs6m1$Drug<- factor(rs6m1$Drug, levels = c('none', 'Methyl', 'Psilo','Drug1','Drug2','Drug'))
        rs6m2$Drug=rep('none', nrow(rs6m2))
        rs6m2$Drug<- factor(rs6m2$Drug, levels = c('none', 'Methyl', 'Psilo','Drug1','Drug2','Drug'))
        # for resting-state 1
	        for (s in 1:length(unique(rs1m2$Subjects))){
                  # this subject
                  subj=unique(rs1m2$Subjects)[s]
            # Subset rows for the current subject
            subjRows <- rs1m2[rs1m2$Subjects == subj, ]
            # Subset rows for the 'Drug' dosage
            subjDrugRows <- subjRows[subjRows$Dosage == 'Drug', ]
              # Identify 'Drug1' and 'Drug2' based on 'OgRow'
              drug1Rows <- which(subjDrugRows$Session == 1)
              drug2Rows <- which(subjDrugRows$Session == 2)
                  subjRows$Drug[drug1Rows] <- 'Drug1'
                  subjRows$Drug[drug2Rows] <- 'Drug2'
                    ### Now pull out methyl vs. psilo
                    DoseCorresp=subjDoseCorresp[subjDoseCorresp$X1==subj,]
                    # if Drug1 = 1 in key, Drug1 = psilo
                    if (subjDrugRows$Drug[1] == 'DrugOnly_CantSeeMe') {
                                subjRows$Drug[subjRows$Drug == 'Drug'] <- 'Psilo'
                      } else if (DoseCorresp$X2 == 1) {
                                  subjRows$Drug[subjRows$Drug == 'Drug1'] <- 'Psilo'
                          subjRows$Drug[subjRows$Drug == 'Drug2'] <- 'Methyl'
                            } else if (DoseCorresp$X2 == 2) {
                                        subjRows$Drug[subjRows$Drug == 'Drug1'] <- 'Methyl'
                              subjRows$Drug[subjRows$Drug == 'Drug2'] <- 'Psilo'
                                }
                      # find where OG row corresponds to methyl
                      MethylRows=subjRows$OgRow[subjRows$Drug=='Methyl']
                      # and psilo
                      PsiloRows=subjRows$OgRow[subjRows$Drug=='Psilo']
                        # set
                        rs1m2$Drug[rs1m2$OgRow==MethylRows[1]]='Methyl'
                        rs1m2$Drug[rs1m2$OgRow==PsiloRows[1]]='Psilo'
        }
        # for resting-state 2
        for (s in 1:length(unique(rs2m2$Subjects))){
                  # this subject
                  subj=unique(rs2m2$Subjects)[s]
            # Subset rows for the current subject
            subjRows <- rs2m2[rs2m2$Subjects == subj, ]
            # Subset rows for the 'Drug' dosage
            subjDrugRows <- subjRows[subjRows$Dosage == 'Drug', ]
              # Identify 'Drug1' and 'Drug2' based on 'OgRow'
              drug1Rows <- which(subjDrugRows$Session == 1)
              drug2Rows <- which(subjDrugRows$Session == 2)
                #  # Assign 'Drug' values to the original rows
                  subjRows$Drug[drug1Rows] <- 'Drug1'
                  subjRows$Drug[drug2Rows] <- 'Drug2'
                    ### Now pull out methyl vs. psilo
                    DoseCorresp=subjDoseCorresp[subjDoseCorresp$X1==subj,]
                    # if Drug1 = 1 in key, Drug1 = psilo
                    if (subjDrugRows$Drug[1] == 'DrugOnly_CantSeeMe') {
                                subjRows$Drug[subjRows$Drug == 'Drug'] <- 'Psilo'
                      } else if (DoseCorresp$X2 == 1) {
                                  subjRows$Drug[subjRows$Drug == 'Drug1'] <- 'Psilo'
                          subjRows$Drug[subjRows$Drug == 'Drug2'] <- 'Methyl'
                            } else if (DoseCorresp$X2 == 2) {
                                        subjRows$Drug[subjRows$Drug == 'Drug1'] <- 'Methyl'
                              subjRows$Drug[subjRows$Drug == 'Drug2'] <- 'Psilo'
                                }
                      # find where OG row corresponds to methyl
                      MethylRows=subjRows$OgRow[subjRows$Drug=='Methyl']
                      # and psilo
                      PsiloRows=subjRows$OgRow[subjRows$Drug=='Psilo']
                        # set
                        rs2m2$Drug[rs2m2$OgRow==MethylRows[1]]='Methyl'
                        rs2m2$Drug[rs2m2$OgRow==PsiloRows[1]]='Psilo'
        }
        # for resting-state 3
        for (s in 1:length(unique(rs3m2$Subjects))){
                  # this subject
                  subj=unique(rs3m2$Subjects)[s]
            # Subset rows for the current subject
            subjRows <- rs3m2[rs3m2$Subjects == subj, ]
            # Subset rows for the 'Drug' dosage
            subjDrugRows <- subjRows[subjRows$Dosage == 'Drug', ]
                # Identify 'Drug1' and 'Drug2' based on 'OgRow'
              drug1Rows <- which(subjDrugRows$Session == 1)
              drug2Rows <- which(subjDrugRows$Session == 2)
                # if min=max, label only 1 drug
                #  # Assign 'Drug' values to the original rows
                  subjRows$Drug[drug1Rows] <- 'Drug1'
                  subjRows$Drug[drug2Rows] <- 'Drug2'
                    #}
                    ### Now pull out methyl vs. psilo
                    DoseCorresp=subjDoseCorresp[subjDoseCorresp$X1==subj,]
                    # if Drug1 = 1 in key, Drug1 = psilo
                    if (subjDrugRows$Drug[1] == 'DrugOnly_cantSeeMe') {
                                subjRows$Drug[subjRows$Drug == 'Drug'] <- 'Psilo'
                      } else if (DoseCorresp$X2 == 1) {
                                  subjRows$Drug[subjRows$Drug == 'Drug1'] <- 'Psilo'
                          subjRows$Drug[subjRows$Drug == 'Drug2'] <- 'Methyl'
                            } else if (DoseCorresp$X2 == 2) {
                                        subjRows$Drug[subjRows$Drug == 'Drug1'] <- 'Methyl'
                              subjRows$Drug[subjRows$Drug == 'Drug2'] <- 'Psilo'
                                }
                      # find where OG row corresponds to methyl
                      MethylRows=subjRows$OgRow[subjRows$Drug=='Methyl']
                      # and psilo
                      PsiloRows=subjRows$OgRow[subjRows$Drug=='Psilo']
                        # set
                        rs3m2$Drug[rs3m2$OgRow==MethylRows[1]]='Methyl'
                        rs3m2$Drug[rs3m2$OgRow==PsiloRows[1]]='Psilo'
        }
        # for resting-state 4
        for (s in 1:length(unique(rs4m2$Subjects))){
                  # this subject
                  subj=unique(rs4m2$Subjects)[s]
            # Subset rows for the current subject
            subjRows <- rs4m2[rs4m2$Subjects == subj, ]
            # Subset rows for the 'Drug' dosage
            subjDrugRows <- subjRows[subjRows$Dosage == 'Drug', ]
              # Identify 'Drug1' and 'Drug2' based on 'OgRow'
              drug1Rows <- which(subjDrugRows$Session == 1)
              drug2Rows <- which(subjDrugRows$Session == 2)
                # if min=max, label only 1 drug
                #  # Assign 'Drug' values to the original rows
                  subjRows$Drug[drug1Rows] <- 'Drug1'
                  subjRows$Drug[drug2Rows] <- 'Drug2'
                    #}
                    ### Now pull out methyl vs. psilo
                    DoseCorresp=subjDoseCorresp[subjDoseCorresp$X1==subj,]
                    # if Drug1 = 1 in key, Drug1 = psilo
                    if (subjDrugRows$Drug[1] == 'DrugOnly_CantSeeMe') {
                                subjRows$Drug[subjRows$Drug == 'Drug'] <- 'Psilo'
                      } else if (DoseCorresp$X2 == 1) {
                                  subjRows$Drug[subjRows$Drug == 'Drug1'] <- 'Psilo'
                          subjRows$Drug[subjRows$Drug == 'Drug2'] <- 'Methyl'
                            } else if (DoseCorresp$X2 == 2) {
                                        subjRows$Drug[subjRows$Drug == 'Drug1'] <- 'Methyl'
                              subjRows$Drug[subjRows$Drug == 'Drug2'] <- 'Psilo'
                                }
                      # find where OG row corresponds to methyl
                      MethylRows=subjRows$OgRow[subjRows$Drug=='Methyl']
                      # and psilo
                      PsiloRows=subjRows$OgRow[subjRows$Drug=='Psilo']
                        # set
                        rs4m2$Drug[rs4m2$OgRow==MethylRows[1]]='Methyl'
                        rs4m2$Drug[rs4m2$OgRow==PsiloRows[1]]='Psilo'
        }
        # for resting-state 5
        for (s in 1:length(unique(rs5m2$Subjects))){
                  # this subject
                  subj=unique(rs5m2$Subjects)[s]
            # Subset rows for the current subject
            subjRows <- rs5m2[rs5m2$Subjects == subj, ]
            # Subset rows for the 'Drug' dosage
            subjDrugRows <- subjRows[subjRows$Dosage == 'Drug', ]
              # Identify 'Drug1' and 'Drug2' based on 'OgRow'
              drug1Rows <- which(subjDrugRows$Session == 1)
              drug2Rows <- which(subjDrugRows$Session == 2)
                #  # Assign 'Drug' values to the original rows
                  subjRows$Drug[drug1Rows] <- 'Drug1'
                  subjRows$Drug[drug2Rows] <- 'Drug2'
                    #}
                    ### Now pull out methyl vs. psilo
                    DoseCorresp=subjDoseCorresp[subjDoseCorresp$X1==subj,]
                    # if Drug1 = 1 in key, Drug1 = psilo
                    if (subjDrugRows$Drug[1] == 'DrugOnly_CantSeeMe') {
                                subjRows$Drug[subjRows$Drug == 'Drug'] <- 'Psilo'
                      } else if (DoseCorresp$X2 == 1) {
                                  subjRows$Drug[subjRows$Drug == 'Drug1'] <- 'Psilo'
                          subjRows$Drug[subjRows$Drug == 'Drug2'] <- 'Methyl'
                            } else if (DoseCorresp$X2 == 2) {
                                        subjRows$Drug[subjRows$Drug == 'Drug1'] <- 'Methyl'
                              subjRows$Drug[subjRows$Drug == 'Drug2'] <- 'Psilo'
                                }
                      # find where OG row corresponds to methyl
                      MethylRows=subjRows$OgRow[subjRows$Drug=='Methyl']
                      # and psilo
                      PsiloRows=subjRows$OgRow[subjRows$Drug=='Psilo']
                        # set
                        rs5m2$Drug[rs5m2$OgRow==MethylRows[1]]='Methyl'
                        rs5m2$Drug[rs5m2$OgRow==PsiloRows[1]]='Psilo'
        }
        # for resting-state 6
        for (s in 1:length(unique(rs6m2$Subjects))){
                  # this subject
                  subj=unique(rs6m2$Subjects)[s]
            # Subset rows for the current subject
            subjRows <- rs6m2[rs6m2$Subjects == subj, ]
            # Subset rows for the 'Drug' dosage
            subjDrugRows <- subjRows[subjRows$Dosage == 'Drug', ]
              # Identify 'Drug1' and 'Drug2' based on 'OgRow'
              drug1Rows <- which(subjDrugRows$Session == 1)
              drug2Rows <- which(subjDrugRows$Session == 2)
                #  # Assign 'Drug' values to the original rows
                  subjRows$Drug[drug1Rows] <- 'Drug1'
                  subjRows$Drug[drug2Rows] <- 'Drug2'
                    ### Now pull out methyl vs. psilo
                    DoseCorresp=subjDoseCorresp[subjDoseCorresp$X1==subj,]
                    # if Drug1 = 1 in key, Drug1 = psilo
                    if (subjDrugRows$Drug[1] == 'DrugOnly_CantSeeMe') {
                                subjRows$Drug[subjRows$Drug == 'Drug'] <- 'Psilo'
                      } else if (DoseCorresp$X2 == 1) {
                                  subjRows$Drug[subjRows$Drug == 'Drug1'] <- 'Psilo'
                          subjRows$Drug[subjRows$Drug == 'Drug2'] <- 'Methyl'
                            } else if (DoseCorresp$X2 == 2) {
                                        subjRows$Drug[subjRows$Drug == 'Drug1'] <- 'Methyl'
                              subjRows$Drug[subjRows$Drug == 'Drug2'] <- 'Psilo'
                                }
                      # find where OG row corresponds to methyl
                      MethylRows=subjRows$OgRow[subjRows$Drug=='Methyl']
                      # and psilo
                      PsiloRows=subjRows$OgRow[subjRows$Drug=='Psilo']
                        # set
                        rs6m2$Drug[rs6m2$OgRow==MethylRows[1]]='Methyl'
                        rs6m2$Drug[rs6m2$OgRow==PsiloRows[1]]='Psilo'
        }

        # combine all
        allScans=rbind(rs1bv,rs1p,rs1m1,rs1m2,rs2bv,rs2p,rs2m1,rs2m2,rs3bv,rs3p,rs3m1,rs3m2,rs4bv,rs4p,rs4m1,rs4m2,rs5bv,rs5p,rs5m1,rs5m2,rs6bv,rs6p,rs6m1,rs6m2)
        # before scans
        Before=rbind(rs1bv,rs2bv,rs3bv,rs4bv,rs5bv,rs6bv)
        Before$Chronology='Before'
        # after
        After=rbind(rs1m1,rs2m1,rs3m1,rs4m1,rs5m1,rs6m1)
        After$Chronology='After'
        B_A=rbind(Before,After)
        # set non-drug conditions to none
        allScans$Drug[allScans$Dosage=='none']='none'
        # retain a legacy subjects naming convention for MEQ use later
        allScans$LegacySubjs=allScans$Subjects
        # final thing is to convert all "before" to "Baseline", "between" to "Between". and "after" to "After"
        allScans$Session <- gsub("before", "Baseline", allScans$Session)
        allScans$Session <- gsub("between", "Between", allScans$Session)
        allScans$Session <- gsub("after", "After", allScans$Session)
        # qc criteria
        allScans=allScans[allScans$RemTRs>250,]
        # model
        allScansNoMeth=allScans[allScans$Drug!='Methyl',]
        fit_lme <- lme(TDProp1 ~ Drug + RemTRs + FD, random = ~ 1 | Subjects, data = allScansNoMeth)
	summaryLME<-summary(fit_lme)
        tdistr_Mag[n]=summaryLME$tTable['DrugPsilo','t-value']
}
print('done with mags')

for (n in 1:numSpins){
        rs1=read.csv(paste0('/scratch/users/apines/data/rs1_Psil_VISTemporalAutoCor_Merged_Spin_',n,'.csv'),header=F)
        rs2=read.csv(paste0('/scratch/users/apines/data/rs2_Psil_VISTemporalAutoCor_Merged_Spin_',n,'.csv'),header=F)
        rs3=read.csv(paste0('/scratch/users/apines/data/rs3_Psil_VISTemporalAutoCor_Merged_Spin_',n,'.csv'),header=F)
        rs4=read.csv(paste0('/scratch/users/apines/data/rs4_Psil_VISTemporalAutoCor_Merged_Spin_',n,'.csv'),header=F)
        rs5=read.csv(paste0('/scratch/users/apines/data/rs5_Psil_VISTemporalAutoCor_Merged_Spin_',n,'.csv'),header=F)
        rs6=read.csv(paste0('/scratch/users/apines/data/rs6_Psil_VISTemporalAutoCor_Merged_Spin_',n,'.csv'),header=F)
        # mimic main analysis code
        rs1$Task='rs'
        rs2$Task='rs2'
        rs3$Task='rs3'
        rs4$Task='rs4'
        rs5$Task='rs5'
        rs6$Task='rs6'
        rs1$OgRow=seq(1:dim(rs1)[1])
        rs2$OgRow=seq(1:dim(rs2)[1])
        rs3$OgRow=seq(1:dim(rs3)[1])
        rs4$OgRow=seq(1:dim(rs4)[1])
        rs5$OgRow=seq(1:dim(rs5)[1])
        rs6$OgRow=seq(1:dim(rs6)[1])
        rs1bv=data.frame(cbind(rs1$V1,rs1$V2,rs1$V3,rs1$V4,rs1$V17,rs1$OgRow,rs1$V21))
        rs1p=data.frame(cbind(rs1$V5,rs1$V6,rs1$V7,rs1$V8,rs1$V18,rs1$OgRow,rs1$V22))
        rs1m1=data.frame(cbind(rs1$V9,rs1$V10,rs1$V11,rs1$V12,rs1$V19,rs1$OgRow,rs1$V23))
        rs1m2=data.frame(cbind(rs1$V13,rs1$V14,rs1$V15,rs1$V16,rs1$V20,rs1$OgRow,rs1$V24))
        colnames(rs1bv)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs','OgRow','FD')
        colnames(rs1p)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs','OgRow','FD')
        colnames(rs1m1)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs','OgRow','FD')
        colnames(rs1m2)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs','OgRow','FD')
        rs2bv=data.frame(cbind(rs2$V1,rs2$V2,rs2$V3,rs2$V4,rs2$V17,rs2$OgRow,rs2$V21))
        rs2p=data.frame(cbind(rs2$V5,rs2$V6,rs2$V7,rs2$V8,rs2$V18,rs2$OgRow,rs2$V22))
        rs2m1=data.frame(cbind(rs2$V9,rs2$V10,rs2$V11,rs2$V12,rs2$V19,rs2$OgRow,rs2$V23))
        rs2m2=data.frame(cbind(rs2$V13,rs2$V14,rs2$V15,rs2$V16,rs2$V20,rs2$OgRow,rs2$V24))
        colnames(rs2bv)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs','OgRow','FD')
        colnames(rs2p)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs','OgRow','FD')
        colnames(rs2m1)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs','OgRow','FD')
        colnames(rs2m2)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs','OgRow','FD')
        rs3bv=data.frame(cbind(rs3$V1,rs3$V2,rs3$V3,rs3$V4,rs3$V17,rs3$OgRow,rs3$V21))
        rs3p=data.frame(cbind(rs3$V5,rs3$V6,rs3$V7,rs3$V8,rs3$V18,rs3$OgRow,rs3$V22))
        rs3m1=data.frame(cbind(rs3$V9,rs3$V10,rs3$V11,rs3$V12,rs3$V19,rs3$OgRow,rs3$V23))
        rs3m2=data.frame(cbind(rs3$V13,rs3$V14,rs3$V15,rs3$V16,rs3$V20,rs3$OgRow,rs3$V24))
        colnames(rs3bv)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs','OgRow','FD')
        colnames(rs3p)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs','OgRow','FD')
        colnames(rs3m1)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs','OgRow','FD')
        colnames(rs3m2)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs','OgRow','FD')
        rs4bv=data.frame(cbind(rs4$V1,rs4$V2,rs4$V3,rs4$V4,rs4$V17,rs4$OgRow,rs4$V21))
        rs4p=data.frame(cbind(rs4$V5,rs4$V6,rs4$V7,rs4$V8,rs4$V18,rs4$OgRow,rs4$V22))
        rs4m1=data.frame(cbind(rs4$V9,rs4$V10,rs4$V11,rs4$V12,rs4$V19,rs4$OgRow,rs4$V23))
        rs4m2=data.frame(cbind(rs4$V13,rs4$V14,rs4$V15,rs4$V16,rs4$V20,rs4$OgRow,rs4$V24))
        colnames(rs4bv)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs','OgRow','FD')
        colnames(rs4p)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs','OgRow','FD')
        colnames(rs4m1)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs','OgRow','FD')
        colnames(rs4m2)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs','OgRow','FD')
        rs5bv=data.frame(cbind(rs5$V1,rs5$V2,rs5$V3,rs5$V4,rs5$V17,rs5$OgRow,rs5$V21))
        rs5p=data.frame(cbind(rs5$V5,rs5$V6,rs5$V7,rs5$V8,rs5$V18,rs5$OgRow,rs5$V22))
        rs5m1=data.frame(cbind(rs5$V9,rs5$V10,rs5$V11,rs5$V12,rs5$V19,rs5$OgRow,rs5$V23))
        rs5m2=data.frame(cbind(rs5$V13,rs5$V14,rs5$V15,rs5$V16,rs5$V20,rs5$OgRow,rs5$V24))
        colnames(rs5bv)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs','OgRow','FD')
        colnames(rs5p)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs','OgRow','FD')
        colnames(rs5m1)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs','OgRow','FD')
        colnames(rs5m2)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs','OgRow','FD')
        rs6bv=data.frame(cbind(rs6$V1,rs6$V2,rs6$V3,rs6$V4,rs6$V17,rs6$OgRow,rs6$V21))
        rs6p=data.frame(cbind(rs6$V5,rs6$V6,rs6$V7,rs6$V8,rs6$V18,rs6$OgRow,rs6$V22))
        rs6m1=data.frame(cbind(rs6$V9,rs6$V10,rs6$V11,rs6$V12,rs6$V19,rs6$OgRow,rs6$V23))
        rs6m2=data.frame(cbind(rs6$V13,rs6$V14,rs6$V15,rs6$V16,rs6$V20,rs6$OgRow,rs6$V24))
        colnames(rs6bv)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs','OgRow','FD')
        colnames(rs6p)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs','OgRow','FD')
        colnames(rs6m1)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs','OgRow','FD')
        colnames(rs6m2)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs','OgRow','FD')
        # get subject IDs subject order csv (should both be identical)
        subjOrder_rs1=read.csv('/oak/stanford/groups/leanew1/users/apines/data/rs1_Psil_propsMerged_subjOrder.csv',blank.lines.skip = FALSE)
        subjOrder_rs2=read.csv('/oak/stanford/groups/leanew1/users/apines/data/rs2_Psil_propsMerged_subjOrder.csv',blank.lines.skip = FALSE)
        subjOrder_rs3=read.csv('/oak/stanford/groups/leanew1/users/apines/data/rs3_Psil_propsMerged_subjOrder.csv',blank.lines.skip = FALSE)
        subjOrder_rs4=read.csv('/oak/stanford/groups/leanew1/users/apines/data/rs4_Psil_propsMerged_subjOrder.csv',blank.lines.skip = FALSE)
        subjOrder_rs5=read.csv('/oak/stanford/groups/leanew1/users/apines/data/rs5_Psil_propsMerged_subjOrder.csv',blank.lines.skip = FALSE)
        subjOrder_rs6=read.csv('/oak/stanford/groups/leanew1/users/apines/data/rs6_Psil_propsMerged_subjOrder.csv',blank.lines.skip = FALSE)
        rs1bv$Subjects=subjOrder_rs1$SubjNameCol_1
        rs2bv$Subjects=subjOrder_rs2$SubjNameCol_1
        rs3bv$Subjects=subjOrder_rs3$SubjNameCol_1
        rs4bv$Subjects=subjOrder_rs4$SubjNameCol_1
        rs5bv$Subjects=subjOrder_rs5$SubjNameCol_1
        rs6bv$Subjects=subjOrder_rs6$SubjNameCol_1
        rs1p$Subjects=subjOrder_rs1$SubjNameCol_1
        rs2p$Subjects=subjOrder_rs2$SubjNameCol_1
        rs3p$Subjects=subjOrder_rs3$SubjNameCol_1
        rs4p$Subjects=subjOrder_rs4$SubjNameCol_1
        rs5p$Subjects=subjOrder_rs5$SubjNameCol_1
        rs6p$Subjects=subjOrder_rs6$SubjNameCol_1
        rs1m1$Subjects=subjOrder_rs1$SubjNameCol_1
        rs2m1$Subjects=subjOrder_rs2$SubjNameCol_1
        rs3m1$Subjects=subjOrder_rs3$SubjNameCol_1
        rs4m1$Subjects=subjOrder_rs4$SubjNameCol_1
        rs5m1$Subjects=subjOrder_rs5$SubjNameCol_1
        rs6m1$Subjects=subjOrder_rs6$SubjNameCol_1
        rs1m2$Subjects=subjOrder_rs1$SubjNameCol_1
        rs2m2$Subjects=subjOrder_rs2$SubjNameCol_1
        rs3m2$Subjects=subjOrder_rs3$SubjNameCol_1
        rs4m2$Subjects=subjOrder_rs4$SubjNameCol_1
        rs5m2$Subjects=subjOrder_rs5$SubjNameCol_1
        rs6m2$Subjects=subjOrder_rs6$SubjNameCol_1
        # add session straight from matlab readout
        subjOrder_rs_ws=read.delim('~/rs_Psil_propsMerged_subjOrder_WithSesh.csv',blank.lines.skip = FALSE,sep=',')
        rs1bv$Session=subjOrder_rs_ws$SubjNameCol_2
        rs2bv$Session=subjOrder_rs_ws$SubjNameCol_2
        rs3bv$Session=subjOrder_rs_ws$SubjNameCol_2
        rs4bv$Session=subjOrder_rs_ws$SubjNameCol_2
        rs5bv$Session=subjOrder_rs_ws$SubjNameCol_2
        rs6bv$Session=subjOrder_rs_ws$SubjNameCol_2
        rs1p$Session=subjOrder_rs_ws$SubjNameCol_2
        rs2p$Session=subjOrder_rs_ws$SubjNameCol_2
        rs3p$Session=subjOrder_rs_ws$SubjNameCol_2
        rs4p$Session=subjOrder_rs_ws$SubjNameCol_2
        rs5p$Session=subjOrder_rs_ws$SubjNameCol_2
        rs6p$Session=subjOrder_rs_ws$SubjNameCol_2
        rs1m1$Session=subjOrder_rs_ws$SubjNameCol_2
        rs2m1$Session=subjOrder_rs_ws$SubjNameCol_2
        rs3m1$Session=subjOrder_rs_ws$SubjNameCol_2
        rs4m1$Session=subjOrder_rs_ws$SubjNameCol_2
        rs5m1$Session=subjOrder_rs_ws$SubjNameCol_2
        rs6m1$Session=subjOrder_rs_ws$SubjNameCol_2
        rs1m2$Session=subjOrder_rs_ws$SubjNameCol_2
        rs2m2$Session=subjOrder_rs_ws$SubjNameCol_2
        rs3m2$Session=subjOrder_rs_ws$SubjNameCol_2
        rs4m2$Session=subjOrder_rs_ws$SubjNameCol_2
        rs5m2$Session=subjOrder_rs_ws$SubjNameCol_2
        rs6m2$Session=subjOrder_rs_ws$SubjNameCol_2
        # add in task (rs to be made equivalent after motion merge)
        rs1bv$Task='rs'
        rs1p$Task='rs'
        rs1m1$Task='rs'
        rs1m2$Task='rs'
        rs2bv$Task='rs2'
        rs2p$Task='rs2'
        rs2m1$Task='rs2'
        rs2m2$Task='rs2'
        rs3bv$Task='rs3'
        rs3p$Task='rs3'
        rs3m1$Task='rs3'
        rs3m2$Task='rs3'
        rs4bv$Task='rs4'
        rs4p$Task='rs4'
        rs4m1$Task='rs4'
        rs4m2$Task='rs4'
        rs5bv$Task='rs5'
        rs5p$Task='rs5'
        rs5m1$Task='rs5'
        rs5m2$Task='rs5'
        rs6bv$Task='rs6'
        rs6p$Task='rs6'
        rs6m1$Task='rs6'
        rs6m2$Task='rs6'
        # add in dosage
        rs1bv$Dosage='none'
        rs1p$Dosage='none'
        rs1m1$Dosage='none'
        rs1m2$Dosage='Drug'
        rs2bv$Dosage='none'
        rs2p$Dosage='none'
        rs2m1$Dosage='none'
        rs2m2$Dosage='Drug'
        rs3bv$Dosage='none'
        rs3p$Dosage='none'
        rs3m1$Dosage='none'
        rs3m2$Dosage='Drug'
        rs4bv$Dosage='none'
        rs4p$Dosage='none'
        rs4m1$Dosage='none'
        rs4m2$Dosage='Drug'
        rs5bv$Dosage='none'
        rs5p$Dosage='none'
        rs5m1$Dosage='none'
        rs5m2$Dosage='Drug'
        rs6bv$Dosage='none'
        rs6p$Dosage='none'
        rs6m1$Dosage='none'
        rs6m2$Dosage='Drug'
        # parse out only existing rows
        rs1bv=rs1bv[rs1bv$TDProp1>0,]
        rs1p=rs1p[rs1p$TDProp1>0,]
        rs1m1=rs1m1[rs1m1$TDProp1>0,]
        rs1m2=rs1m2[rs1m2$TDProp1>0,]
        rs2bv=rs2bv[rs2bv$TDProp1>0,]
        rs2p=rs2p[rs2p$TDProp1>0,]
        rs2m1=rs2m1[rs2m1$TDProp1>0,]
        rs2m2=rs2m2[rs2m2$TDProp1>0,]
        rs3bv=rs3bv[rs3bv$TDProp1>0,]
        rs3p=rs3p[rs3p$TDProp1>0,]
        rs3m1=rs3m1[rs3m1$TDProp1>0,]
        rs3m2=rs3m2[rs3m2$TDProp1>0,]
        rs4bv=rs4bv[rs4bv$TDProp1>0,]
        rs4p=rs4p[rs4p$TDProp1>0,]
        rs4m1=rs4m1[rs4m1$TDProp1>0,]
        rs4m2=rs4m2[rs4m2$TDProp1>0,]
        rs5bv=rs5bv[rs5bv$TDProp1>0,]
        rs5p=rs5p[rs5p$TDProp1>0,]
        rs5m1=rs5m1[rs5m1$TDProp1>0,]
        rs5m2=rs5m2[rs5m2$TDProp1>0,]
        rs6bv=rs6bv[rs6bv$TDProp1>0,]
        rs6p=rs6p[rs6p$TDProp1>0,]
        rs6m1=rs6m1[rs6m1$TDProp1>0,]
        rs6m2=rs6m2[rs6m2$TDProp1>0,]
        # decode Methylphenidate vs. psilocybin for each PT
        subjDoseCorresp=read.csv('~/subjSeshDoseCorresp_psilo.csv',header=F)
        subjDoseCorresp=data.frame(t(subjDoseCorresp))
        subjDoseCorresp$X2<-as.numeric(subjDoseCorresp$X2)
        # initialize new drug column in all dfs
        rs1bv$Drug <- rep('none', nrow(rs1bv))
        rs1bv$Drug <- factor(rs1bv$Drug, levels = c('none', 'Methyl', 'Psilo','Drug1','Drug2','Drug'))
        rs2bv$Drug=rep('none', nrow(rs2bv))
        rs2bv$Drug<- factor(rs2bv$Drug, levels = c('none', 'Methyl', 'Psilo','Drug1','Drug2','Drug'))
        rs1p$Drug=rep('none', nrow(rs1p))
        rs1p$Drug<- factor(rs1p$Drug, levels = c('none', 'Methyl', 'Psilo','Drug1','Drug2','Drug'))
        rs2p$Drug=rep('none', nrow(rs2p))
        rs2p$Drug<- factor(rs2p$Drug, levels = c('none', 'Methyl', 'Psilo','Drug1','Drug2','Drug'))
        rs1m1$Drug=rep('none', nrow(rs1m1))
        rs1m1$Drug<- factor(rs1m1$Drug, levels = c('none', 'Methyl', 'Psilo','Drug1','Drug2','Drug'))
        rs2m1$Drug=rep('none', nrow(rs2m1))
        rs2m1$Drug<- factor(rs2m1$Drug, levels = c('none', 'Methyl', 'Psilo','Drug1','Drug2','Drug'))
        rs1m2$Drug=rep('none', nrow(rs1m2))
        rs1m2$Drug<- factor(rs1m2$Drug, levels = c('none', 'Methyl', 'Psilo','Drug1','Drug2','Drug'))
        rs2m2$Drug=rep('none', nrow(rs2m2))
        rs2m2$Drug<- factor(rs2m2$Drug, levels = c('none', 'Methyl', 'Psilo','Drug1','Drug2','Drug'))
        rs3bv$Drug=rep('none', nrow(rs3bv))
        rs3bv$Drug<- factor(rs3bv$Drug, levels = c('none', 'Methyl', 'Psilo','Drug1','Drug2','Drug'))
        rs3p$Drug=rep('none', nrow(rs3p))
        rs3p$Drug<- factor(rs3p$Drug, levels = c('none', 'Methyl', 'Psilo','Drug1','Drug2','Drug'))
        rs3m1$Drug=rep('none', nrow(rs3m1))
        rs3m1$Drug<- factor(rs3m1$Drug, levels = c('none', 'Methyl', 'Psilo','Drug1','Drug2','Drug'))
        rs3m2$Drug=rep('none', nrow(rs3m2))
        rs3m2$Drug<- factor(rs3m2$Drug, levels = c('none', 'Methyl', 'Psilo','Drug1','Drug2','Drug'))
        rs4bv$Drug=rep('none', nrow(rs4bv))
        rs4bv$Drug<- factor(rs4bv$Drug, levels = c('none', 'Methyl', 'Psilo','Drug1','Drug2','Drug'))
        rs4p$Drug=rep('none', nrow(rs4p))
        rs4p$Drug<- factor(rs4p$Drug, levels = c('none', 'Methyl', 'Psilo','Drug1','Drug2','Drug'))
        rs4m1$Drug=rep('none', nrow(rs4m1))
        rs4m1$Drug<- factor(rs4m1$Drug, levels = c('none', 'Methyl', 'Psilo','Drug1','Drug2','Drug'))
        rs4m2$Drug=rep('none', nrow(rs4m2))
        rs4m2$Drug<- factor(rs4m2$Drug, levels = c('none', 'Methyl', 'Psilo','Drug1','Drug2','Drug'))
        rs5bv$Drug=rep('none', nrow(rs5bv))
        rs5bv$Drug<- factor(rs5bv$Drug, levels = c('none', 'Methyl', 'Psilo','Drug1','Drug2','Drug'))
        rs5p$Drug=rep('none', nrow(rs5p))
        rs5p$Drug<- factor(rs5p$Drug, levels = c('none', 'Methyl', 'Psilo','Drug1','Drug2','Drug'))
        rs5m1$Drug=rep('none', nrow(rs5m1))
        rs5m1$Drug<- factor(rs5m1$Drug, levels = c('none', 'Methyl', 'Psilo','Drug1','Drug2','Drug'))
        rs5m2$Drug=rep('none', nrow(rs5m2))
        rs5m2$Drug<- factor(rs5m2$Drug, levels = c('none', 'Methyl', 'Psilo','Drug1','Drug2','Drug'))
        rs6bv$Drug=rep('none', nrow(rs6bv))
        rs6bv$Drug<- factor(rs6bv$Drug, levels = c('none', 'Methyl', 'Psilo','Drug1','Drug2','Drug'))
        rs6p$Drug=rep('none', nrow(rs6p))
        rs6p$Drug<- factor(rs6p$Drug, levels = c('none', 'Methyl', 'Psilo','Drug1','Drug2','Drug'))
        rs6m1$Drug=rep('none', nrow(rs6m1))
        rs6m1$Drug<- factor(rs6m1$Drug, levels = c('none', 'Methyl', 'Psilo','Drug1','Drug2','Drug'))
        rs6m2$Drug=rep('none', nrow(rs6m2))
        rs6m2$Drug<- factor(rs6m2$Drug, levels = c('none', 'Methyl', 'Psilo','Drug1','Drug2','Drug'))
        # for resting-state 1
        for (s in 1:length(unique(rs1m2$Subjects))){
                  # this subject
                  subj=unique(rs1m2$Subjects)[s]
            # Subset rows for the current subject
            subjRows <- rs1m2[rs1m2$Subjects == subj, ]
            # Subset rows for the 'Drug' dosage
            subjDrugRows <- subjRows[subjRows$Dosage == 'Drug', ]
              # Identify 'Drug1' and 'Drug2' based on 'OgRow'
              drug1Rows <- which(subjDrugRows$Session == 1)
              drug2Rows <- which(subjDrugRows$Session == 2)
                  subjRows$Drug[drug1Rows] <- 'Drug1'
                  subjRows$Drug[drug2Rows] <- 'Drug2'
                    ### Now pull out methyl vs. psilo
                    DoseCorresp=subjDoseCorresp[subjDoseCorresp$X1==subj,]
                    # if Drug1 = 1 in key, Drug1 = psilo
                    if (subjDrugRows$Drug[1] == 'DrugOnly_CantSeeMe') {
                                subjRows$Drug[subjRows$Drug == 'Drug'] <- 'Psilo'
                      } else if (DoseCorresp$X2 == 1) {
                                  subjRows$Drug[subjRows$Drug == 'Drug1'] <- 'Psilo'
                          subjRows$Drug[subjRows$Drug == 'Drug2'] <- 'Methyl'
                            } else if (DoseCorresp$X2 == 2) {
                                        subjRows$Drug[subjRows$Drug == 'Drug1'] <- 'Methyl'
                              subjRows$Drug[subjRows$Drug == 'Drug2'] <- 'Psilo'
                                }
                      # find where OG row corresponds to methyl
                      MethylRows=subjRows$OgRow[subjRows$Drug=='Methyl']
                      # and psilo
                      PsiloRows=subjRows$OgRow[subjRows$Drug=='Psilo']
                        # set
                        rs1m2$Drug[rs1m2$OgRow==MethylRows[1]]='Methyl'
                        rs1m2$Drug[rs1m2$OgRow==PsiloRows[1]]='Psilo'
        }
        # for resting-state 2
        for (s in 1:length(unique(rs2m2$Subjects))){
                  # this subject
                  subj=unique(rs2m2$Subjects)[s]
            # Subset rows for the current subject
            subjRows <- rs2m2[rs2m2$Subjects == subj, ]
            # Subset rows for the 'Drug' dosage
            subjDrugRows <- subjRows[subjRows$Dosage == 'Drug', ]
              # Identify 'Drug1' and 'Drug2' based on 'OgRow'
              drug1Rows <- which(subjDrugRows$Session == 1)
              drug2Rows <- which(subjDrugRows$Session == 2)
                #  # Assign 'Drug' values to the original rows
                  subjRows$Drug[drug1Rows] <- 'Drug1'
                  subjRows$Drug[drug2Rows] <- 'Drug2'
                    ### Now pull out methyl vs. psilo
                    DoseCorresp=subjDoseCorresp[subjDoseCorresp$X1==subj,]
                    # if Drug1 = 1 in key, Drug1 = psilo
                    if (subjDrugRows$Drug[1] == 'DrugOnly_CantSeeMe') {
                                subjRows$Drug[subjRows$Drug == 'Drug'] <- 'Psilo'
                      } else if (DoseCorresp$X2 == 1) {
                                  subjRows$Drug[subjRows$Drug == 'Drug1'] <- 'Psilo'
                          subjRows$Drug[subjRows$Drug == 'Drug2'] <- 'Methyl'
                            } else if (DoseCorresp$X2 == 2) {
                                        subjRows$Drug[subjRows$Drug == 'Drug1'] <- 'Methyl'
                              subjRows$Drug[subjRows$Drug == 'Drug2'] <- 'Psilo'
                                }
                      # find where OG row corresponds to methyl
                      MethylRows=subjRows$OgRow[subjRows$Drug=='Methyl']
                      # and psilo
                      PsiloRows=subjRows$OgRow[subjRows$Drug=='Psilo']
                        # set
                        rs2m2$Drug[rs2m2$OgRow==MethylRows[1]]='Methyl'
                        rs2m2$Drug[rs2m2$OgRow==PsiloRows[1]]='Psilo'
        }
        # for resting-state 3
        for (s in 1:length(unique(rs3m2$Subjects))){
                  # this subject
                  subj=unique(rs3m2$Subjects)[s]
            # Subset rows for the current subject
            subjRows <- rs3m2[rs3m2$Subjects == subj, ]
            # Subset rows for the 'Drug' dosage
            subjDrugRows <- subjRows[subjRows$Dosage == 'Drug', ]
                # Identify 'Drug1' and 'Drug2' based on 'OgRow'
              drug1Rows <- which(subjDrugRows$Session == 1)
              drug2Rows <- which(subjDrugRows$Session == 2)
                # if min=max, label only 1 drug
                #  # Assign 'Drug' values to the original rows
                  subjRows$Drug[drug1Rows] <- 'Drug1'
                  subjRows$Drug[drug2Rows] <- 'Drug2'
                    #}
                    ### Now pull out methyl vs. psilo
                    DoseCorresp=subjDoseCorresp[subjDoseCorresp$X1==subj,]
                    # if Drug1 = 1 in key, Drug1 = psilo
                    if (subjDrugRows$Drug[1] == 'DrugOnly_cantSeeMe') {
                                subjRows$Drug[subjRows$Drug == 'Drug'] <- 'Psilo'
                      } else if (DoseCorresp$X2 == 1) {
                                  subjRows$Drug[subjRows$Drug == 'Drug1'] <- 'Psilo'
                          subjRows$Drug[subjRows$Drug == 'Drug2'] <- 'Methyl'
                            } else if (DoseCorresp$X2 == 2) {
                                        subjRows$Drug[subjRows$Drug == 'Drug1'] <- 'Methyl'
                              subjRows$Drug[subjRows$Drug == 'Drug2'] <- 'Psilo'
                                }
                      # find where OG row corresponds to methyl
                      MethylRows=subjRows$OgRow[subjRows$Drug=='Methyl']
                      # and psilo
                      PsiloRows=subjRows$OgRow[subjRows$Drug=='Psilo']
                        # set
                        rs3m2$Drug[rs3m2$OgRow==MethylRows[1]]='Methyl'
                        rs3m2$Drug[rs3m2$OgRow==PsiloRows[1]]='Psilo'
        }
        # for resting-state 4
        for (s in 1:length(unique(rs4m2$Subjects))){
                  # this subject
                  subj=unique(rs4m2$Subjects)[s]
            # Subset rows for the current subject
            subjRows <- rs4m2[rs4m2$Subjects == subj, ]
            # Subset rows for the 'Drug' dosage
            subjDrugRows <- subjRows[subjRows$Dosage == 'Drug', ]
              # Identify 'Drug1' and 'Drug2' based on 'OgRow'
              drug1Rows <- which(subjDrugRows$Session == 1)
              drug2Rows <- which(subjDrugRows$Session == 2)
                # if min=max, label only 1 drug
                #  # Assign 'Drug' values to the original rows
                  subjRows$Drug[drug1Rows] <- 'Drug1'
                  subjRows$Drug[drug2Rows] <- 'Drug2'
                    #}
                    ### Now pull out methyl vs. psilo
                    DoseCorresp=subjDoseCorresp[subjDoseCorresp$X1==subj,]
                    # if Drug1 = 1 in key, Drug1 = psilo
                    if (subjDrugRows$Drug[1] == 'DrugOnly_CantSeeMe') {
                                subjRows$Drug[subjRows$Drug == 'Drug'] <- 'Psilo'
                      } else if (DoseCorresp$X2 == 1) {
                                  subjRows$Drug[subjRows$Drug == 'Drug1'] <- 'Psilo'
                          subjRows$Drug[subjRows$Drug == 'Drug2'] <- 'Methyl'
                            } else if (DoseCorresp$X2 == 2) {
                                        subjRows$Drug[subjRows$Drug == 'Drug1'] <- 'Methyl'
                              subjRows$Drug[subjRows$Drug == 'Drug2'] <- 'Psilo'
                                }
                      # find where OG row corresponds to methyl
                      MethylRows=subjRows$OgRow[subjRows$Drug=='Methyl']
                      # and psilo
                      PsiloRows=subjRows$OgRow[subjRows$Drug=='Psilo']
                        # set
                        rs4m2$Drug[rs4m2$OgRow==MethylRows[1]]='Methyl'
                        rs4m2$Drug[rs4m2$OgRow==PsiloRows[1]]='Psilo'
        }
        # for resting-state 5
        for (s in 1:length(unique(rs5m2$Subjects))){
                  # this subject
                  subj=unique(rs5m2$Subjects)[s]
            # Subset rows for the current subject
            subjRows <- rs5m2[rs5m2$Subjects == subj, ]
            # Subset rows for the 'Drug' dosage
            subjDrugRows <- subjRows[subjRows$Dosage == 'Drug', ]
              # Identify 'Drug1' and 'Drug2' based on 'OgRow'
              drug1Rows <- which(subjDrugRows$Session == 1)
              drug2Rows <- which(subjDrugRows$Session == 2)
                #  # Assign 'Drug' values to the original rows
                  subjRows$Drug[drug1Rows] <- 'Drug1'
                  subjRows$Drug[drug2Rows] <- 'Drug2'
                    #}
                    ### Now pull out methyl vs. psilo
                    DoseCorresp=subjDoseCorresp[subjDoseCorresp$X1==subj,]
                    # if Drug1 = 1 in key, Drug1 = psilo
                    if (subjDrugRows$Drug[1] == 'DrugOnly_CantSeeMe') {
                                subjRows$Drug[subjRows$Drug == 'Drug'] <- 'Psilo'
                      } else if (DoseCorresp$X2 == 1) {
                                  subjRows$Drug[subjRows$Drug == 'Drug1'] <- 'Psilo'
                          subjRows$Drug[subjRows$Drug == 'Drug2'] <- 'Methyl'
                            } else if (DoseCorresp$X2 == 2) {
                                        subjRows$Drug[subjRows$Drug == 'Drug1'] <- 'Methyl'
                              subjRows$Drug[subjRows$Drug == 'Drug2'] <- 'Psilo'
                                }
                      # find where OG row corresponds to methyl
                      MethylRows=subjRows$OgRow[subjRows$Drug=='Methyl']
                      # and psilo
                      PsiloRows=subjRows$OgRow[subjRows$Drug=='Psilo']
                        # set
                        rs5m2$Drug[rs5m2$OgRow==MethylRows[1]]='Methyl'
                        rs5m2$Drug[rs5m2$OgRow==PsiloRows[1]]='Psilo'
        }
        # for resting-state 6
        for (s in 1:length(unique(rs6m2$Subjects))){
                  # this subject
                  subj=unique(rs6m2$Subjects)[s]
            # Subset rows for the current subject
            subjRows <- rs6m2[rs6m2$Subjects == subj, ]
            # Subset rows for the 'Drug' dosage
            subjDrugRows <- subjRows[subjRows$Dosage == 'Drug', ]
              # Identify 'Drug1' and 'Drug2' based on 'OgRow'
              drug1Rows <- which(subjDrugRows$Session == 1)
              drug2Rows <- which(subjDrugRows$Session == 2)
                #  # Assign 'Drug' values to the original rows
                  subjRows$Drug[drug1Rows] <- 'Drug1'
                  subjRows$Drug[drug2Rows] <- 'Drug2'
                    ### Now pull out methyl vs. psilo
                    DoseCorresp=subjDoseCorresp[subjDoseCorresp$X1==subj,]
                    # if Drug1 = 1 in key, Drug1 = psilo
                    if (subjDrugRows$Drug[1] == 'DrugOnly_CantSeeMe') {
                                subjRows$Drug[subjRows$Drug == 'Drug'] <- 'Psilo'
                      } else if (DoseCorresp$X2 == 1) {
                                  subjRows$Drug[subjRows$Drug == 'Drug1'] <- 'Psilo'
                          subjRows$Drug[subjRows$Drug == 'Drug2'] <- 'Methyl'
                            } else if (DoseCorresp$X2 == 2) {
                                        subjRows$Drug[subjRows$Drug == 'Drug1'] <- 'Methyl'
                              subjRows$Drug[subjRows$Drug == 'Drug2'] <- 'Psilo'
                                }
                      # find where OG row corresponds to methyl
                      MethylRows=subjRows$OgRow[subjRows$Drug=='Methyl']
                      # and psilo
                      PsiloRows=subjRows$OgRow[subjRows$Drug=='Psilo']
                        # set
                        rs6m2$Drug[rs6m2$OgRow==MethylRows[1]]='Methyl'
                        rs6m2$Drug[rs6m2$OgRow==PsiloRows[1]]='Psilo'
        }

        # combine all
        allScans=rbind(rs1bv,rs1p,rs1m1,rs1m2,rs2bv,rs2p,rs2m1,rs2m2,rs3bv,rs3p,rs3m1,rs3m2,rs4bv,rs4p,rs4m1,rs4m2,rs5bv,rs5p,rs5m1,rs5m2,rs6bv,rs6p,rs6m1,rs6m2)
        # before scans
        Before=rbind(rs1bv,rs2bv,rs3bv,rs4bv,rs5bv,rs6bv)
        Before$Chronology='Before'
        # after
        After=rbind(rs1m1,rs2m1,rs3m1,rs4m1,rs5m1,rs6m1)
        After$Chronology='After'
        B_A=rbind(Before,After)
        # set non-drug conditions to none
        allScans$Drug[allScans$Dosage=='none']='none'
        # retain a legacy subjects naming convention for MEQ use later
        allScans$LegacySubjs=allScans$Subjects
        # final thing is to convert all "before" to "Baseline", "between" to "Between". and "after" to "After"
        allScans$Session <- gsub("before", "Baseline", allScans$Session)
        allScans$Session <- gsub("between", "Between", allScans$Session)
        allScans$Session <- gsub("after", "After", allScans$Session)
        # qc criteria
        allScans=allScans[allScans$RemTRs>250,]
        # model
        allScansNoMeth=allScans[allScans$Drug!='Methyl',]
        fit_lme <- lme(TDProp1 ~ Drug + RemTRs + FD, random = ~ 1 | Subjects, data = allScansNoMeth)
 	summaryLME<-summary(fit_lme)
 	tdistr_TA[n]=summaryLME$tTable['DrugPsilo','t-value']
}
print('done with DMN T autocor')
######################
# dmn fc
################
for (n in 1:numSpins){
        rs1=read.csv(paste0('/scratch/users/apines/data/rs1_Psil_VISSegMerged_Spin_',n,'.csv'),header=F)
        rs2=read.csv(paste0('/scratch/users/apines/data/rs2_Psil_VISSegMerged_Spin_',n,'.csv'),header=F)
        rs3=read.csv(paste0('/scratch/users/apines/data/rs3_Psil_VISSegMerged_Spin_',n,'.csv'),header=F)
        rs4=read.csv(paste0('/scratch/users/apines/data/rs4_Psil_VISSegMerged_Spin_',n,'.csv'),header=F)
        rs5=read.csv(paste0('/scratch/users/apines/data/rs5_Psil_VISSegMerged_Spin_',n,'.csv'),header=F)
        rs6=read.csv(paste0('/scratch/users/apines/data/rs6_Psil_VISSegMerged_Spin_',n,'.csv'),header=F)
        # mimic main analysis code
        rs1$Task='rs'
        rs2$Task='rs2'
        rs3$Task='rs3'
        rs4$Task='rs4'
        rs5$Task='rs5'
        rs6$Task='rs6'
        rs1$OgRow=seq(1:dim(rs1)[1])
        rs2$OgRow=seq(1:dim(rs2)[1])
        rs3$OgRow=seq(1:dim(rs3)[1])
        rs4$OgRow=seq(1:dim(rs4)[1])
        rs5$OgRow=seq(1:dim(rs5)[1])
        rs6$OgRow=seq(1:dim(rs6)[1])
        rs1bv=data.frame(cbind(rs1$V1,rs1$V2,rs1$V3,rs1$V4,rs1$V17,rs1$OgRow,rs1$V21))
        rs1p=data.frame(cbind(rs1$V5,rs1$V6,rs1$V7,rs1$V8,rs1$V18,rs1$OgRow,rs1$V22))
        rs1m1=data.frame(cbind(rs1$V9,rs1$V10,rs1$V11,rs1$V12,rs1$V19,rs1$OgRow,rs1$V23))
        rs1m2=data.frame(cbind(rs1$V13,rs1$V14,rs1$V15,rs1$V16,rs1$V20,rs1$OgRow,rs1$V24))
        colnames(rs1bv)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs','OgRow','FD')
        colnames(rs1p)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs','OgRow','FD')
        colnames(rs1m1)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs','OgRow','FD')
        colnames(rs1m2)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs','OgRow','FD')
        rs2bv=data.frame(cbind(rs2$V1,rs2$V2,rs2$V3,rs2$V4,rs2$V17,rs2$OgRow,rs2$V21))
        rs2p=data.frame(cbind(rs2$V5,rs2$V6,rs2$V7,rs2$V8,rs2$V18,rs2$OgRow,rs2$V22))
        rs2m1=data.frame(cbind(rs2$V9,rs2$V10,rs2$V11,rs2$V12,rs2$V19,rs2$OgRow,rs2$V23))
        rs2m2=data.frame(cbind(rs2$V13,rs2$V14,rs2$V15,rs2$V16,rs2$V20,rs2$OgRow,rs2$V24))
        colnames(rs2bv)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs','OgRow','FD')
        colnames(rs2p)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs','OgRow','FD')
        colnames(rs2m1)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs','OgRow','FD')
        colnames(rs2m2)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs','OgRow','FD')
        rs3bv=data.frame(cbind(rs3$V1,rs3$V2,rs3$V3,rs3$V4,rs3$V17,rs3$OgRow,rs3$V21))
        rs3p=data.frame(cbind(rs3$V5,rs3$V6,rs3$V7,rs3$V8,rs3$V18,rs3$OgRow,rs3$V22))
        rs3m1=data.frame(cbind(rs3$V9,rs3$V10,rs3$V11,rs3$V12,rs3$V19,rs3$OgRow,rs3$V23))
        rs3m2=data.frame(cbind(rs3$V13,rs3$V14,rs3$V15,rs3$V16,rs3$V20,rs3$OgRow,rs3$V24))
        colnames(rs3bv)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs','OgRow','FD')
        colnames(rs3p)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs','OgRow','FD')
        colnames(rs3m1)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs','OgRow','FD')
        colnames(rs3m2)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs','OgRow','FD')
        rs4bv=data.frame(cbind(rs4$V1,rs4$V2,rs4$V3,rs4$V4,rs4$V17,rs4$OgRow,rs4$V21))
        rs4p=data.frame(cbind(rs4$V5,rs4$V6,rs4$V7,rs4$V8,rs4$V18,rs4$OgRow,rs4$V22))
        rs4m1=data.frame(cbind(rs4$V9,rs4$V10,rs4$V11,rs4$V12,rs4$V19,rs4$OgRow,rs4$V23))
        rs4m2=data.frame(cbind(rs4$V13,rs4$V14,rs4$V15,rs4$V16,rs4$V20,rs4$OgRow,rs4$V24))
        colnames(rs4bv)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs','OgRow','FD')
        colnames(rs4p)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs','OgRow','FD')
        colnames(rs4m1)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs','OgRow','FD')
        colnames(rs4m2)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs','OgRow','FD')
        rs5bv=data.frame(cbind(rs5$V1,rs5$V2,rs5$V3,rs5$V4,rs5$V17,rs5$OgRow,rs5$V21))
        rs5p=data.frame(cbind(rs5$V5,rs5$V6,rs5$V7,rs5$V8,rs5$V18,rs5$OgRow,rs5$V22))
        rs5m1=data.frame(cbind(rs5$V9,rs5$V10,rs5$V11,rs5$V12,rs5$V19,rs5$OgRow,rs5$V23))
        rs5m2=data.frame(cbind(rs5$V13,rs5$V14,rs5$V15,rs5$V16,rs5$V20,rs5$OgRow,rs5$V24))
        colnames(rs5bv)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs','OgRow','FD')
        colnames(rs5p)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs','OgRow','FD')
        colnames(rs5m1)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs','OgRow','FD')
        colnames(rs5m2)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs','OgRow','FD')
        rs6bv=data.frame(cbind(rs6$V1,rs6$V2,rs6$V3,rs6$V4,rs6$V17,rs6$OgRow,rs6$V21))
        rs6p=data.frame(cbind(rs6$V5,rs6$V6,rs6$V7,rs6$V8,rs6$V18,rs6$OgRow,rs6$V22))
        rs6m1=data.frame(cbind(rs6$V9,rs6$V10,rs6$V11,rs6$V12,rs6$V19,rs6$OgRow,rs6$V23))
        rs6m2=data.frame(cbind(rs6$V13,rs6$V14,rs6$V15,rs6$V16,rs6$V20,rs6$OgRow,rs6$V24))
        colnames(rs6bv)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs','OgRow','FD')
        colnames(rs6p)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs','OgRow','FD')
        colnames(rs6m1)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs','OgRow','FD')
        colnames(rs6m2)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs','OgRow','FD')
        # get subject IDs subject order csv (should both be identical)
        subjOrder_rs1=read.csv('/oak/stanford/groups/leanew1/users/apines/data/rs1_Psil_propsMerged_subjOrder.csv',blank.lines.skip = FALSE)
        subjOrder_rs2=read.csv('/oak/stanford/groups/leanew1/users/apines/data/rs2_Psil_propsMerged_subjOrder.csv',blank.lines.skip = FALSE)
        subjOrder_rs3=read.csv('/oak/stanford/groups/leanew1/users/apines/data/rs3_Psil_propsMerged_subjOrder.csv',blank.lines.skip = FALSE)
        subjOrder_rs4=read.csv('/oak/stanford/groups/leanew1/users/apines/data/rs4_Psil_propsMerged_subjOrder.csv',blank.lines.skip = FALSE)
        subjOrder_rs5=read.csv('/oak/stanford/groups/leanew1/users/apines/data/rs5_Psil_propsMerged_subjOrder.csv',blank.lines.skip = FALSE)
        subjOrder_rs6=read.csv('/oak/stanford/groups/leanew1/users/apines/data/rs6_Psil_propsMerged_subjOrder.csv',blank.lines.skip = FALSE)
        rs1bv$Subjects=subjOrder_rs1$SubjNameCol_1
        rs2bv$Subjects=subjOrder_rs2$SubjNameCol_1
        rs3bv$Subjects=subjOrder_rs3$SubjNameCol_1
        rs4bv$Subjects=subjOrder_rs4$SubjNameCol_1
        rs5bv$Subjects=subjOrder_rs5$SubjNameCol_1
        rs6bv$Subjects=subjOrder_rs6$SubjNameCol_1
        rs1p$Subjects=subjOrder_rs1$SubjNameCol_1
        rs2p$Subjects=subjOrder_rs2$SubjNameCol_1
        rs3p$Subjects=subjOrder_rs3$SubjNameCol_1
        rs4p$Subjects=subjOrder_rs4$SubjNameCol_1
        rs5p$Subjects=subjOrder_rs5$SubjNameCol_1
        rs6p$Subjects=subjOrder_rs6$SubjNameCol_1
        rs1m1$Subjects=subjOrder_rs1$SubjNameCol_1
        rs2m1$Subjects=subjOrder_rs2$SubjNameCol_1
        rs3m1$Subjects=subjOrder_rs3$SubjNameCol_1
        rs4m1$Subjects=subjOrder_rs4$SubjNameCol_1
        rs5m1$Subjects=subjOrder_rs5$SubjNameCol_1
        rs6m1$Subjects=subjOrder_rs6$SubjNameCol_1
        rs1m2$Subjects=subjOrder_rs1$SubjNameCol_1
        rs2m2$Subjects=subjOrder_rs2$SubjNameCol_1
        rs3m2$Subjects=subjOrder_rs3$SubjNameCol_1
        rs4m2$Subjects=subjOrder_rs4$SubjNameCol_1
        rs5m2$Subjects=subjOrder_rs5$SubjNameCol_1
        rs6m2$Subjects=subjOrder_rs6$SubjNameCol_1
        # add session straight from matlab readout
        subjOrder_rs_ws=read.delim('~/rs_Psil_propsMerged_subjOrder_WithSesh.csv',blank.lines.skip = FALSE,sep=',')
        rs1bv$Session=subjOrder_rs_ws$SubjNameCol_2
        rs2bv$Session=subjOrder_rs_ws$SubjNameCol_2
        rs3bv$Session=subjOrder_rs_ws$SubjNameCol_2
        rs4bv$Session=subjOrder_rs_ws$SubjNameCol_2
        rs5bv$Session=subjOrder_rs_ws$SubjNameCol_2
        rs6bv$Session=subjOrder_rs_ws$SubjNameCol_2
        rs1p$Session=subjOrder_rs_ws$SubjNameCol_2
        rs2p$Session=subjOrder_rs_ws$SubjNameCol_2
        rs3p$Session=subjOrder_rs_ws$SubjNameCol_2
        rs4p$Session=subjOrder_rs_ws$SubjNameCol_2
        rs5p$Session=subjOrder_rs_ws$SubjNameCol_2
        rs6p$Session=subjOrder_rs_ws$SubjNameCol_2
        rs1m1$Session=subjOrder_rs_ws$SubjNameCol_2
        rs2m1$Session=subjOrder_rs_ws$SubjNameCol_2
        rs3m1$Session=subjOrder_rs_ws$SubjNameCol_2
        rs4m1$Session=subjOrder_rs_ws$SubjNameCol_2
        rs5m1$Session=subjOrder_rs_ws$SubjNameCol_2
        rs6m1$Session=subjOrder_rs_ws$SubjNameCol_2
        rs1m2$Session=subjOrder_rs_ws$SubjNameCol_2
        rs2m2$Session=subjOrder_rs_ws$SubjNameCol_2
        rs3m2$Session=subjOrder_rs_ws$SubjNameCol_2
        rs4m2$Session=subjOrder_rs_ws$SubjNameCol_2
        rs5m2$Session=subjOrder_rs_ws$SubjNameCol_2
        rs6m2$Session=subjOrder_rs_ws$SubjNameCol_2
        # add in task (rs to be made equivalent after motion merge)
        rs1bv$Task='rs'
        rs1p$Task='rs'
        rs1m1$Task='rs'
        rs1m2$Task='rs'
        rs2bv$Task='rs2'
        rs2p$Task='rs2'
        rs2m1$Task='rs2'
        rs2m2$Task='rs2'
        rs3bv$Task='rs3'
        rs3p$Task='rs3'
        rs3m1$Task='rs3'
        rs3m2$Task='rs3'
        rs4bv$Task='rs4'
        rs4p$Task='rs4'
        rs4m1$Task='rs4'
        rs4m2$Task='rs4'
        rs5bv$Task='rs5'
        rs5p$Task='rs5'
        rs5m1$Task='rs5'
        rs5m2$Task='rs5'
        rs6bv$Task='rs6'
        rs6p$Task='rs6'
        rs6m1$Task='rs6'
        rs6m2$Task='rs6'
        # add in dosage
        rs1bv$Dosage='none'
        rs1p$Dosage='none'
        rs1m1$Dosage='none'
        rs1m2$Dosage='Drug'
        rs2bv$Dosage='none'
        rs2p$Dosage='none'
        rs2m1$Dosage='none'
        rs2m2$Dosage='Drug'
        rs3bv$Dosage='none'
        rs3p$Dosage='none'
        rs3m1$Dosage='none'
        rs3m2$Dosage='Drug'
        rs4bv$Dosage='none'
        rs4p$Dosage='none'
        rs4m1$Dosage='none'
        rs4m2$Dosage='Drug'
        rs5bv$Dosage='none'
        rs5p$Dosage='none'
        rs5m1$Dosage='none'
        rs5m2$Dosage='Drug'
        rs6bv$Dosage='none'
        rs6p$Dosage='none'
        rs6m1$Dosage='none'
        rs6m2$Dosage='Drug'
        # parse out only existing rows
        rs1bv=rs1bv[rs1bv$TDProp1>0,]
        rs1p=rs1p[rs1p$TDProp1>0,]
        rs1m1=rs1m1[rs1m1$TDProp1>0,]
        rs1m2=rs1m2[rs1m2$TDProp1>0,]
        rs2bv=rs2bv[rs2bv$TDProp1>0,]
        rs2p=rs2p[rs2p$TDProp1>0,]
        rs2m1=rs2m1[rs2m1$TDProp1>0,]
        rs2m2=rs2m2[rs2m2$TDProp1>0,]
        rs3bv=rs3bv[rs3bv$TDProp1>0,]
        rs3p=rs3p[rs3p$TDProp1>0,]
        rs3m1=rs3m1[rs3m1$TDProp1>0,]
        rs3m2=rs3m2[rs3m2$TDProp1>0,]
        rs4bv=rs4bv[rs4bv$TDProp1>0,]
        rs4p=rs4p[rs4p$TDProp1>0,]
        rs4m1=rs4m1[rs4m1$TDProp1>0,]
        rs4m2=rs4m2[rs4m2$TDProp1>0,]
        rs5bv=rs5bv[rs5bv$TDProp1>0,]
        rs5p=rs5p[rs5p$TDProp1>0,]
        rs5m1=rs5m1[rs5m1$TDProp1>0,]
        rs5m2=rs5m2[rs5m2$TDProp1>0,]
        rs6bv=rs6bv[rs6bv$TDProp1>0,]
        rs6p=rs6p[rs6p$TDProp1>0,]
        rs6m1=rs6m1[rs6m1$TDProp1>0,]
        rs6m2=rs6m2[rs6m2$TDProp1>0,]
        # decode Methylphenidate vs. psilocybin for each PT
        subjDoseCorresp=read.csv('~/subjSeshDoseCorresp_psilo.csv',header=F)
        subjDoseCorresp=data.frame(t(subjDoseCorresp))
        subjDoseCorresp$X2<-as.numeric(subjDoseCorresp$X2)
        # initialize new drug column in all dfs
        rs1bv$Drug <- rep('none', nrow(rs1bv))
        rs1bv$Drug <- factor(rs1bv$Drug, levels = c('none', 'Methyl', 'Psilo','Drug1','Drug2','Drug'))
        rs2bv$Drug=rep('none', nrow(rs2bv))
        rs2bv$Drug<- factor(rs2bv$Drug, levels = c('none', 'Methyl', 'Psilo','Drug1','Drug2','Drug'))
        rs1p$Drug=rep('none', nrow(rs1p))
        rs1p$Drug<- factor(rs1p$Drug, levels = c('none', 'Methyl', 'Psilo','Drug1','Drug2','Drug'))
        rs2p$Drug=rep('none', nrow(rs2p))
        rs2p$Drug<- factor(rs2p$Drug, levels = c('none', 'Methyl', 'Psilo','Drug1','Drug2','Drug'))
        rs1m1$Drug=rep('none', nrow(rs1m1))
        rs1m1$Drug<- factor(rs1m1$Drug, levels = c('none', 'Methyl', 'Psilo','Drug1','Drug2','Drug'))
        rs2m1$Drug=rep('none', nrow(rs2m1))
        rs2m1$Drug<- factor(rs2m1$Drug, levels = c('none', 'Methyl', 'Psilo','Drug1','Drug2','Drug'))
        rs1m2$Drug=rep('none', nrow(rs1m2))
        rs1m2$Drug<- factor(rs1m2$Drug, levels = c('none', 'Methyl', 'Psilo','Drug1','Drug2','Drug'))
        rs2m2$Drug=rep('none', nrow(rs2m2))
        rs2m2$Drug<- factor(rs2m2$Drug, levels = c('none', 'Methyl', 'Psilo','Drug1','Drug2','Drug'))
        rs3bv$Drug=rep('none', nrow(rs3bv))
        rs3bv$Drug<- factor(rs3bv$Drug, levels = c('none', 'Methyl', 'Psilo','Drug1','Drug2','Drug'))
        rs3p$Drug=rep('none', nrow(rs3p))
        rs3p$Drug<- factor(rs3p$Drug, levels = c('none', 'Methyl', 'Psilo','Drug1','Drug2','Drug'))
        rs3m1$Drug=rep('none', nrow(rs3m1))
        rs3m1$Drug<- factor(rs3m1$Drug, levels = c('none', 'Methyl', 'Psilo','Drug1','Drug2','Drug'))
        rs3m2$Drug=rep('none', nrow(rs3m2))
        rs3m2$Drug<- factor(rs3m2$Drug, levels = c('none', 'Methyl', 'Psilo','Drug1','Drug2','Drug'))
        rs4bv$Drug=rep('none', nrow(rs4bv))
        rs4bv$Drug<- factor(rs4bv$Drug, levels = c('none', 'Methyl', 'Psilo','Drug1','Drug2','Drug'))
        rs4p$Drug=rep('none', nrow(rs4p))
        rs4p$Drug<- factor(rs4p$Drug, levels = c('none', 'Methyl', 'Psilo','Drug1','Drug2','Drug'))
        rs4m1$Drug=rep('none', nrow(rs4m1))
        rs4m1$Drug<- factor(rs4m1$Drug, levels = c('none', 'Methyl', 'Psilo','Drug1','Drug2','Drug'))
        rs4m2$Drug=rep('none', nrow(rs4m2))
        rs4m2$Drug<- factor(rs4m2$Drug, levels = c('none', 'Methyl', 'Psilo','Drug1','Drug2','Drug'))
        rs5bv$Drug=rep('none', nrow(rs5bv))
        rs5bv$Drug<- factor(rs5bv$Drug, levels = c('none', 'Methyl', 'Psilo','Drug1','Drug2','Drug'))
        rs5p$Drug=rep('none', nrow(rs5p))
        rs5p$Drug<- factor(rs5p$Drug, levels = c('none', 'Methyl', 'Psilo','Drug1','Drug2','Drug'))
        rs5m1$Drug=rep('none', nrow(rs5m1))
        rs5m1$Drug<- factor(rs5m1$Drug, levels = c('none', 'Methyl', 'Psilo','Drug1','Drug2','Drug'))
        rs5m2$Drug=rep('none', nrow(rs5m2))
        rs5m2$Drug<- factor(rs5m2$Drug, levels = c('none', 'Methyl', 'Psilo','Drug1','Drug2','Drug'))
        rs6bv$Drug=rep('none', nrow(rs6bv))
        rs6bv$Drug<- factor(rs6bv$Drug, levels = c('none', 'Methyl', 'Psilo','Drug1','Drug2','Drug'))
        rs6p$Drug=rep('none', nrow(rs6p))
        rs6p$Drug<- factor(rs6p$Drug, levels = c('none', 'Methyl', 'Psilo','Drug1','Drug2','Drug'))
        rs6m1$Drug=rep('none', nrow(rs6m1))
        rs6m1$Drug<- factor(rs6m1$Drug, levels = c('none', 'Methyl', 'Psilo','Drug1','Drug2','Drug'))
        rs6m2$Drug=rep('none', nrow(rs6m2))
        rs6m2$Drug<- factor(rs6m2$Drug, levels = c('none', 'Methyl', 'Psilo','Drug1','Drug2','Drug'))
        # for resting-state 1
                for (s in 1:length(unique(rs1m2$Subjects))){
                  # this subject
                  subj=unique(rs1m2$Subjects)[s]
            # Subset rows for the current subject
            subjRows <- rs1m2[rs1m2$Subjects == subj, ]
            # Subset rows for the 'Drug' dosage
            subjDrugRows <- subjRows[subjRows$Dosage == 'Drug', ]
              # Identify 'Drug1' and 'Drug2' based on 'OgRow'
              drug1Rows <- which(subjDrugRows$Session == 1)
              drug2Rows <- which(subjDrugRows$Session == 2)
                  subjRows$Drug[drug1Rows] <- 'Drug1'
                  subjRows$Drug[drug2Rows] <- 'Drug2'
                    ### Now pull out methyl vs. psilo
                    DoseCorresp=subjDoseCorresp[subjDoseCorresp$X1==subj,]
                    # if Drug1 = 1 in key, Drug1 = psilo
                    if (subjDrugRows$Drug[1] == 'DrugOnly_CantSeeMe') {
                                subjRows$Drug[subjRows$Drug == 'Drug'] <- 'Psilo'
                      } else if (DoseCorresp$X2 == 1) {
                                  subjRows$Drug[subjRows$Drug == 'Drug1'] <- 'Psilo'
                          subjRows$Drug[subjRows$Drug == 'Drug2'] <- 'Methyl'
                            } else if (DoseCorresp$X2 == 2) {
                                        subjRows$Drug[subjRows$Drug == 'Drug1'] <- 'Methyl'
                              subjRows$Drug[subjRows$Drug == 'Drug2'] <- 'Psilo'
                                }
                      # find where OG row corresponds to methyl
                      MethylRows=subjRows$OgRow[subjRows$Drug=='Methyl']
                      # and psilo
                      PsiloRows=subjRows$OgRow[subjRows$Drug=='Psilo']
                        # set
                        rs1m2$Drug[rs1m2$OgRow==MethylRows[1]]='Methyl'
                        rs1m2$Drug[rs1m2$OgRow==PsiloRows[1]]='Psilo'
        }
        # for resting-state 2
        for (s in 1:length(unique(rs2m2$Subjects))){
                  # this subject
                  subj=unique(rs2m2$Subjects)[s]
            # Subset rows for the current subject
            subjRows <- rs2m2[rs2m2$Subjects == subj, ]
            # Subset rows for the 'Drug' dosage
            subjDrugRows <- subjRows[subjRows$Dosage == 'Drug', ]
              # Identify 'Drug1' and 'Drug2' based on 'OgRow'
              drug1Rows <- which(subjDrugRows$Session == 1)
              drug2Rows <- which(subjDrugRows$Session == 2)
                #  # Assign 'Drug' values to the original rows
                  subjRows$Drug[drug1Rows] <- 'Drug1'
                  subjRows$Drug[drug2Rows] <- 'Drug2'
                    ### Now pull out methyl vs. psilo
                    DoseCorresp=subjDoseCorresp[subjDoseCorresp$X1==subj,]
                    # if Drug1 = 1 in key, Drug1 = psilo
                    if (subjDrugRows$Drug[1] == 'DrugOnly_CantSeeMe') {
                                subjRows$Drug[subjRows$Drug == 'Drug'] <- 'Psilo'
                      } else if (DoseCorresp$X2 == 1) {
                                  subjRows$Drug[subjRows$Drug == 'Drug1'] <- 'Psilo'
                          subjRows$Drug[subjRows$Drug == 'Drug2'] <- 'Methyl'
                            } else if (DoseCorresp$X2 == 2) {
                                        subjRows$Drug[subjRows$Drug == 'Drug1'] <- 'Methyl'
                              subjRows$Drug[subjRows$Drug == 'Drug2'] <- 'Psilo'
                                }
                      # find where OG row corresponds to methyl
                      MethylRows=subjRows$OgRow[subjRows$Drug=='Methyl']
                      # and psilo
                      PsiloRows=subjRows$OgRow[subjRows$Drug=='Psilo']
                        # set
                        rs2m2$Drug[rs2m2$OgRow==MethylRows[1]]='Methyl'
                        rs2m2$Drug[rs2m2$OgRow==PsiloRows[1]]='Psilo'
        }
        # for resting-state 3
        for (s in 1:length(unique(rs3m2$Subjects))){
                  # this subject
                  subj=unique(rs3m2$Subjects)[s]
            # Subset rows for the current subject
            subjRows <- rs3m2[rs3m2$Subjects == subj, ]
            # Subset rows for the 'Drug' dosage
            subjDrugRows <- subjRows[subjRows$Dosage == 'Drug', ]
                # Identify 'Drug1' and 'Drug2' based on 'OgRow'
              drug1Rows <- which(subjDrugRows$Session == 1)
              drug2Rows <- which(subjDrugRows$Session == 2)
                # if min=max, label only 1 drug
                #  # Assign 'Drug' values to the original rows
                  subjRows$Drug[drug1Rows] <- 'Drug1'
                  subjRows$Drug[drug2Rows] <- 'Drug2'
                    #}
                    ### Now pull out methyl vs. psilo
                    DoseCorresp=subjDoseCorresp[subjDoseCorresp$X1==subj,]
                    # if Drug1 = 1 in key, Drug1 = psilo
                    if (subjDrugRows$Drug[1] == 'DrugOnly_cantSeeMe') {
                                subjRows$Drug[subjRows$Drug == 'Drug'] <- 'Psilo'
                      } else if (DoseCorresp$X2 == 1) {
                                  subjRows$Drug[subjRows$Drug == 'Drug1'] <- 'Psilo'
                          subjRows$Drug[subjRows$Drug == 'Drug2'] <- 'Methyl'
                            } else if (DoseCorresp$X2 == 2) {
                                        subjRows$Drug[subjRows$Drug == 'Drug1'] <- 'Methyl'
                              subjRows$Drug[subjRows$Drug == 'Drug2'] <- 'Psilo'
                                }
                      # find where OG row corresponds to methyl
                      MethylRows=subjRows$OgRow[subjRows$Drug=='Methyl']
                      # and psilo
                      PsiloRows=subjRows$OgRow[subjRows$Drug=='Psilo']
                        # set
                        rs3m2$Drug[rs3m2$OgRow==MethylRows[1]]='Methyl'
                        rs3m2$Drug[rs3m2$OgRow==PsiloRows[1]]='Psilo'
        }
        # for resting-state 4
        for (s in 1:length(unique(rs4m2$Subjects))){
                  # this subject
                  subj=unique(rs4m2$Subjects)[s]
          print(subj)
            # Subset rows for the current subject
            subjRows <- rs4m2[rs4m2$Subjects == subj, ]
            # Subset rows for the 'Drug' dosage
            subjDrugRows <- subjRows[subjRows$Dosage == 'Drug', ]
              # Identify 'Drug1' and 'Drug2' based on 'OgRow'
              drug1Rows <- which(subjDrugRows$Session == 1)
              drug2Rows <- which(subjDrugRows$Session == 2)
                # if min=max, label only 1 drug
                #  # Assign 'Drug' values to the original rows
                  subjRows$Drug[drug1Rows] <- 'Drug1'
                  subjRows$Drug[drug2Rows] <- 'Drug2'
                    #}
                    ### Now pull out methyl vs. psilo
                    DoseCorresp=subjDoseCorresp[subjDoseCorresp$X1==subj,]
                    # if Drug1 = 1 in key, Drug1 = psilo
                    if (subjDrugRows$Drug[1] == 'DrugOnly_CantSeeMe') {
                                subjRows$Drug[subjRows$Drug == 'Drug'] <- 'Psilo'
                      } else if (DoseCorresp$X2 == 1) {
                                  subjRows$Drug[subjRows$Drug == 'Drug1'] <- 'Psilo'
                          subjRows$Drug[subjRows$Drug == 'Drug2'] <- 'Methyl'
                            } else if (DoseCorresp$X2 == 2) {
                                        subjRows$Drug[subjRows$Drug == 'Drug1'] <- 'Methyl'
                              subjRows$Drug[subjRows$Drug == 'Drug2'] <- 'Psilo'
                                }
                      # find where OG row corresponds to methyl
                      MethylRows=subjRows$OgRow[subjRows$Drug=='Methyl']
                      # and psilo
                      PsiloRows=subjRows$OgRow[subjRows$Drug=='Psilo']
                        # set
                        rs4m2$Drug[rs4m2$OgRow==MethylRows[1]]='Methyl'
                        rs4m2$Drug[rs4m2$OgRow==PsiloRows[1]]='Psilo'
        }
        # for resting-state 5
        for (s in 1:length(unique(rs5m2$Subjects))){
                  # this subject
                  subj=unique(rs5m2$Subjects)[s]
            # Subset rows for the current subject
            subjRows <- rs5m2[rs5m2$Subjects == subj, ]
            # Subset rows for the 'Drug' dosage
            subjDrugRows <- subjRows[subjRows$Dosage == 'Drug', ]
              # Identify 'Drug1' and 'Drug2' based on 'OgRow'
              drug1Rows <- which(subjDrugRows$Session == 1)
              drug2Rows <- which(subjDrugRows$Session == 2)
                #  # Assign 'Drug' values to the original rows
                  subjRows$Drug[drug1Rows] <- 'Drug1'
                  subjRows$Drug[drug2Rows] <- 'Drug2'
                    #}
                    ### Now pull out methyl vs. psilo
                    DoseCorresp=subjDoseCorresp[subjDoseCorresp$X1==subj,]
                    # if Drug1 = 1 in key, Drug1 = psilo
                    if (subjDrugRows$Drug[1] == 'DrugOnly_CantSeeMe') {
                                subjRows$Drug[subjRows$Drug == 'Drug'] <- 'Psilo'
                      } else if (DoseCorresp$X2 == 1) {
                                  subjRows$Drug[subjRows$Drug == 'Drug1'] <- 'Psilo'
                          subjRows$Drug[subjRows$Drug == 'Drug2'] <- 'Methyl'
                            } else if (DoseCorresp$X2 == 2) {
                                        subjRows$Drug[subjRows$Drug == 'Drug1'] <- 'Methyl'
                              subjRows$Drug[subjRows$Drug == 'Drug2'] <- 'Psilo'
                                }
                      # find where OG row corresponds to methyl
                      MethylRows=subjRows$OgRow[subjRows$Drug=='Methyl']
                      # and psilo
                      PsiloRows=subjRows$OgRow[subjRows$Drug=='Psilo']
                        # set
                        rs5m2$Drug[rs5m2$OgRow==MethylRows[1]]='Methyl'
                        rs5m2$Drug[rs5m2$OgRow==PsiloRows[1]]='Psilo'
        }
        # for resting-state 6
        for (s in 1:length(unique(rs6m2$Subjects))){
                  # this subject
                  subj=unique(rs6m2$Subjects)[s]
            # Subset rows for the current subject
            subjRows <- rs6m2[rs6m2$Subjects == subj, ]
            # Subset rows for the 'Drug' dosage
            subjDrugRows <- subjRows[subjRows$Dosage == 'Drug', ]
              # Identify 'Drug1' and 'Drug2' based on 'OgRow'
              drug1Rows <- which(subjDrugRows$Session == 1)
              drug2Rows <- which(subjDrugRows$Session == 2)
                #  # Assign 'Drug' values to the original rows
                  subjRows$Drug[drug1Rows] <- 'Drug1'
                  subjRows$Drug[drug2Rows] <- 'Drug2'
                    ### Now pull out methyl vs. psilo
                    DoseCorresp=subjDoseCorresp[subjDoseCorresp$X1==subj,]
                    # if Drug1 = 1 in key, Drug1 = psilo
                    if (subjDrugRows$Drug[1] == 'DrugOnly_CantSeeMe') {
                                subjRows$Drug[subjRows$Drug == 'Drug'] <- 'Psilo'
                      } else if (DoseCorresp$X2 == 1) {
                                  subjRows$Drug[subjRows$Drug == 'Drug1'] <- 'Psilo'
                          subjRows$Drug[subjRows$Drug == 'Drug2'] <- 'Methyl'
                            } else if (DoseCorresp$X2 == 2) {
                                        subjRows$Drug[subjRows$Drug == 'Drug1'] <- 'Methyl'
                              subjRows$Drug[subjRows$Drug == 'Drug2'] <- 'Psilo'
                                }
                      # find where OG row corresponds to methyl
                      MethylRows=subjRows$OgRow[subjRows$Drug=='Methyl']
                      # and psilo
                      PsiloRows=subjRows$OgRow[subjRows$Drug=='Psilo']
                        # set
                        rs6m2$Drug[rs6m2$OgRow==MethylRows[1]]='Methyl'
                        rs6m2$Drug[rs6m2$OgRow==PsiloRows[1]]='Psilo'
        }

        # combine all
        allScans=rbind(rs1bv,rs1p,rs1m1,rs1m2,rs2bv,rs2p,rs2m1,rs2m2,rs3bv,rs3p,rs3m1,rs3m2,rs4bv,rs4p,rs4m1,rs4m2,rs5bv,rs5p,rs5m1,rs5m2,rs6bv,rs6p,rs6m1,rs6m2)
        # before scans
        Before=rbind(rs1bv,rs2bv,rs3bv,rs4bv,rs5bv,rs6bv)
        Before$Chronology='Before'
        # after
        After=rbind(rs1m1,rs2m1,rs3m1,rs4m1,rs5m1,rs6m1)
        After$Chronology='After'
        B_A=rbind(Before,After)
        # set non-drug conditions to none
        allScans$Drug[allScans$Dosage=='none']='none'
        # retain a legacy subjects naming convention for MEQ use later
        allScans$LegacySubjs=allScans$Subjects
        # final thing is to convert all "before" to "Baseline", "between" to "Between". and "after" to "After"
        allScans$Session <- gsub("before", "Baseline", allScans$Session)
        allScans$Session <- gsub("between", "Between", allScans$Session)
        allScans$Session <- gsub("after", "After", allScans$Session)
        # qc criteria
        allScans=allScans[allScans$RemTRs>250,]
        # model
        allScansNoMeth=allScans[allScans$Drug!='Methyl',]
        fit_lme <- lme(TDProp1 ~ Drug + RemTRs + FD, random = ~ 1 | Subjects, data = allScansNoMeth)
        summaryLME<-summary(fit_lme)
	tdistr_Seg[n]=summaryLME$tTable['DrugPsilo','t-value']
}
print('done with DMN seg')
outspins=data.frame(tdistr_BUP,tdistr_Mag,tdistr_Seg,tdistr_TA)
saveRDS(outspins,'~/SpunTDistributions_Psil_VIS.rds')
