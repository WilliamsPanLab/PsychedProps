Psilocybin DMN FC
================
2024-02-03

``` r
library(reshape2)
library(ggplot2)
library(visreg)
library(nlme)
```

``` r
# FC measurements. Note that variable names are the SAME as in the main psil.rmd for equivalence, but refer to DMN segregation in this markdown.
rs1=read.csv('~/Downloads/rs1_Psil_DMNSegMerged.csv',header=F)
rs2=read.csv('~/Downloads/rs2_Psil_DMNSegMerged.csv',header=F)
rs3=read.csv('~/Downloads/rs3_Psil_DMNSegMerged.csv',header=F)
rs4=read.csv('~/Downloads/rs4_Psil_DMNSegMerged.csv',header=F)
rs5=read.csv('~/Downloads/rs5_Psil_DMNSegMerged.csv',header=F)
rs6=read.csv('~/Downloads/rs6_Psil_DMNSegMerged.csv',header=F)
# set colnames
#colnames(rs1)=c('bvProp','pProp','m1Prop','m2Prop','bvTRs','pTRs','m1TRs','m2TRs')
#colnames(rs2)=c('bvProp','pProp','m1Prop','m2Prop','bvTRs','pTRs','m1TRs','m2TRs')
#colnames(emo)=c('bvProp','pProp','m1Prop','m2Prop','bvTRs','pTRs','m1TRs','m2TRs')
#colnames(gambling)=c('bvProp','pProp','m1Prop','m2Prop','bvTRs','pTRs','m1TRs','m2TRs')
#colnames(wm)=c('bvProp','pProp','m1Prop','m2Prop','bvTRs','pTRs','m1TRs','m2TRs')
rs1$Task='rs'
rs2$Task='rs2'
rs3$Task='rs3'
rs4$Task='rs4'
rs5$Task='rs5'
rs6$Task='rs6'

# add an OG row column to determine which came first/second chronologically
rs1$OgRow=seq(1:dim(rs1)[1])
rs2$OgRow=seq(1:dim(rs2)[1])
rs3$OgRow=seq(1:dim(rs3)[1])
rs4$OgRow=seq(1:dim(rs4)[1])
rs5$OgRow=seq(1:dim(rs5)[1])
rs6$OgRow=seq(1:dim(rs6)[1])

# this is going to be ugly but simple
# NOTE VARIABLE NAMES ARE FOR EQUIVALENCE WITH MDMA PROC: 
# BV = BASELINE
# P = BETWEEN
# M1 = AFTER
# M2 = DRUG
# manually pair columns as sep. observations of baseline, placebo, 80, 120mg
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
subjOrder_rs1=read.delim('~/rs1_Psil_propsMerged_subjOrder.csv',blank.lines.skip = FALSE)
subjOrder_rs2=read.delim('~/rs2_Psil_propsMerged_subjOrder.csv',blank.lines.skip = FALSE)
subjOrder_rs3=read.delim('~/rs3_Psil_propsMerged_subjOrder.csv',blank.lines.skip = FALSE)
subjOrder_rs4=read.delim('~/rs4_Psil_propsMerged_subjOrder.csv',blank.lines.skip = FALSE)
subjOrder_rs5=read.delim('~/rs5_Psil_propsMerged_subjOrder.csv',blank.lines.skip = FALSE)
subjOrder_rs6=read.delim('~/rs6_Psil_propsMerged_subjOrder.csv',blank.lines.skip = FALSE)
rs1bv$Subjects=subjOrder_rs1$SubjNameCol
rs2bv$Subjects=subjOrder_rs2$SubjNameCol
rs3bv$Subjects=subjOrder_rs3$SubjNameCol
rs4bv$Subjects=subjOrder_rs4$SubjNameCol
rs5bv$Subjects=subjOrder_rs5$SubjNameCol
rs6bv$Subjects=subjOrder_rs6$SubjNameCol
rs1p$Subjects=subjOrder_rs1$SubjNameCol
rs2p$Subjects=subjOrder_rs2$SubjNameCol
rs3p$Subjects=subjOrder_rs3$SubjNameCol
rs4p$Subjects=subjOrder_rs4$SubjNameCol
rs5p$Subjects=subjOrder_rs5$SubjNameCol
rs6p$Subjects=subjOrder_rs6$SubjNameCol
rs1m1$Subjects=subjOrder_rs1$SubjNameCol
rs2m1$Subjects=subjOrder_rs2$SubjNameCol
rs3m1$Subjects=subjOrder_rs3$SubjNameCol
rs4m1$Subjects=subjOrder_rs4$SubjNameCol
rs5m1$Subjects=subjOrder_rs5$SubjNameCol
rs6m1$Subjects=subjOrder_rs6$SubjNameCol
rs1m2$Subjects=subjOrder_rs1$SubjNameCol
rs2m2$Subjects=subjOrder_rs2$SubjNameCol
rs3m2$Subjects=subjOrder_rs3$SubjNameCol
rs4m2$Subjects=subjOrder_rs4$SubjNameCol
rs5m2$Subjects=subjOrder_rs5$SubjNameCol
rs6m2$Subjects=subjOrder_rs6$SubjNameCol
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
```

``` r
# change to account for psy_pfm structure
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
```

``` r
# BEFORE/AFTER VERSION
# change to account for psy_pfm structure
# add in dosage
rs1bv$Dosage='before'
rs1p$Dosage='between'
rs1m1$Dosage='after'
rs1m2$Dosage='Drug'
#
rs2bv$Dosage='before'
rs2p$Dosage='between'
rs2m1$Dosage='after'
rs2m2$Dosage='Drug'
#
rs3bv$Dosage='before'
rs3p$Dosage='between'
rs3m1$Dosage='after'
rs3m2$Dosage='Drug'
#
rs4bv$Dosage='before'
rs4p$Dosage='between'
rs4m1$Dosage='after'
rs4m2$Dosage='Drug'
#
rs5bv$Dosage='before'
rs5p$Dosage='between'
rs5m1$Dosage='after'
rs5m2$Dosage='Drug'
#
rs6bv$Dosage='before'
rs6p$Dosage='between'
rs6m1$Dosage='after'
rs6m2$Dosage='Drug'
```

``` r
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
```

``` r
# decode Methylphenidate vs. psilocybin for each PT
subjDoseCorresp=read.csv('~/subjSeshDoseCorresp_psilo.csv',header=F)
```

    ## Warning in read.table(file = file, header = header, sep = sep, quote = quote, :
    ## incomplete final line found by readTableHeader on
    ## '~/subjSeshDoseCorresp_psilo.csv'

``` r
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
  print(subj)
  # Subset rows for the current subject
  subjRows <- rs1m2[rs1m2$Subjects == subj, ]
  # Subset rows for the 'Drug' dosage
  subjDrugRows <- subjRows[subjRows$Dosage == 'Drug', ]
  # Identify 'Drug1' and 'Drug2' based on 'OgRow'
  drug1Rows <- subjDrugRows[subjDrugRows$OgRow == min(subjDrugRows$OgRow), ]
  drug2Rows <- subjDrugRows[subjDrugRows$OgRow == max(subjDrugRows$OgRow), ]
  # if min=max, label only 1 drug
  if (min(subjDrugRows$OgRow)==max(subjDrugRows$OgRow)){
    subjDrugRows$Drug='DrugOnly'
    subjRows$Drug[subjDrugRows$OgRow==subjDrugRows$OgRow]='Drug'
  } else {
    # Assign 'Drug' values to the original rows
    subjRows$Drug[subjRows$OgRow %in% drug1Rows$OgRow] <- 'Drug1'
    subjRows$Drug[subjRows$OgRow %in% drug2Rows$OgRow] <- 'Drug2'
  }
  ### Now pull out methyl vs. psilo
  DoseCorresp=subjDoseCorresp[subjDoseCorresp$X1==subj,]
  # if Drug1 = 1 in key, Drug1 = psilo
  if (subjDrugRows$Drug[1] == 'DrugOnly') {
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
```

    ## [1] "PS03"
    ## [1] "PS16"
    ## [1] "PS18"
    ## [1] "PS19"
    ## [1] "PS24"
    ## [1] "PS93"
    ## [1] "PS96"
    ## [1] "PS98"
    ## [1] "PS21"

``` r
# for resting-state 2
for (s in 1:length(unique(rs2m2$Subjects))){
  # this subject
  subj=unique(rs2m2$Subjects)[s]
  print(subj)
  # Subset rows for the current subject
  subjRows <- rs2m2[rs2m2$Subjects == subj, ]
  # Subset rows for the 'Drug' dosage
  subjDrugRows <- subjRows[subjRows$Dosage == 'Drug', ]
  # Identify 'Drug1' and 'Drug2' based on 'OgRow'
  drug1Rows <- subjDrugRows[subjDrugRows$OgRow == min(subjDrugRows$OgRow), ]
  drug2Rows <- subjDrugRows[subjDrugRows$OgRow == max(subjDrugRows$OgRow), ]
  # if min=max, label only 1 drug
  if (min(subjDrugRows$OgRow)==max(subjDrugRows$OgRow)){
    subjDrugRows$Drug='DrugOnly'
    subjRows$Drug[subjDrugRows$OgRow==subjDrugRows$OgRow]='Drug'
  } else {
    # Assign 'Drug' values to the original rows
    subjRows$Drug[subjRows$OgRow %in% drug1Rows$OgRow] <- 'Drug1'
    subjRows$Drug[subjRows$OgRow %in% drug2Rows$OgRow] <- 'Drug2'
  }
  ### Now pull out methyl vs. psilo
  DoseCorresp=subjDoseCorresp[subjDoseCorresp$X1==subj,]
  # if Drug1 = 1 in key, Drug1 = psilo
  if (subjDrugRows$Drug[1] == 'DrugOnly') {
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
```

    ## [1] "PS03"
    ## [1] "PS18"
    ## [1] "PS19"
    ## [1] "PS24"
    ## [1] "PS93"
    ## [1] "PS96"
    ## [1] "PS98"
    ## [1] "PS16"
    ## [1] "PS21"

``` r
# for resting-state 3
for (s in 1:length(unique(rs3m2$Subjects))){
  # this subject
  subj=unique(rs3m2$Subjects)[s]
  print(subj)
  # Subset rows for the current subject
  subjRows <- rs3m2[rs3m2$Subjects == subj, ]
  # Subset rows for the 'Drug' dosage
  subjDrugRows <- subjRows[subjRows$Dosage == 'Drug', ]
  # Identify 'Drug1' and 'Drug2' based on 'OgRow'
  drug1Rows <- subjDrugRows[subjDrugRows$OgRow == min(subjDrugRows$OgRow), ]
  drug2Rows <- subjDrugRows[subjDrugRows$OgRow == max(subjDrugRows$OgRow), ]
  # if min=max, label only 1 drug
  if (min(subjDrugRows$OgRow)==max(subjDrugRows$OgRow)){
    subjDrugRows$Drug='DrugOnly'
    subjRows$Drug[subjDrugRows$OgRow==subjDrugRows$OgRow]='Drug'
  } else {
    # Assign 'Drug' values to the original rows
    subjRows$Drug[subjRows$OgRow %in% drug1Rows$OgRow] <- 'Drug1'
    subjRows$Drug[subjRows$OgRow %in% drug2Rows$OgRow] <- 'Drug2'
  }
  ### Now pull out methyl vs. psilo
  DoseCorresp=subjDoseCorresp[subjDoseCorresp$X1==subj,]
  # if Drug1 = 1 in key, Drug1 = psilo
  if (subjDrugRows$Drug[1] == 'DrugOnly') {
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
```

    ## [1] "PS18"
    ## [1] "PS19"
    ## [1] "PS21"
    ## [1] "PS24"
    ## [1] "PS93"
    ## [1] "PS96"

``` r
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
  drug1Rows <- subjDrugRows[subjDrugRows$OgRow == min(subjDrugRows$OgRow), ]
  drug2Rows <- subjDrugRows[subjDrugRows$OgRow == max(subjDrugRows$OgRow), ]
  # if min=max, label only 1 drug
  if (min(subjDrugRows$OgRow)==max(subjDrugRows$OgRow)){
    subjDrugRows$Drug='DrugOnly'
    subjRows$Drug[subjDrugRows$OgRow==subjDrugRows$OgRow]='Drug'
  } else {
    # Assign 'Drug' values to the original rows
    subjRows$Drug[subjRows$OgRow %in% drug1Rows$OgRow] <- 'Drug1'
    subjRows$Drug[subjRows$OgRow %in% drug2Rows$OgRow] <- 'Drug2'
  }
  ### Now pull out methyl vs. psilo
  DoseCorresp=subjDoseCorresp[subjDoseCorresp$X1==subj,]
  # if Drug1 = 1 in key, Drug1 = psilo
  if (subjDrugRows$Drug[1] == 'DrugOnly') {
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
```

    ## [1] "PS19"
    ## [1] "PS24"
    ## [1] "PS93"

``` r
# for resting-state 5
for (s in 1:length(unique(rs5m2$Subjects))){
  # this subject
  subj=unique(rs5m2$Subjects)[s]
  print(subj)
  # Subset rows for the current subject
  subjRows <- rs5m2[rs5m2$Subjects == subj, ]
  # Subset rows for the 'Drug' dosage
  subjDrugRows <- subjRows[subjRows$Dosage == 'Drug', ]
  # Identify 'Drug1' and 'Drug2' based on 'OgRow'
  drug1Rows <- subjDrugRows[subjDrugRows$OgRow == min(subjDrugRows$OgRow), ]
  drug2Rows <- subjDrugRows[subjDrugRows$OgRow == max(subjDrugRows$OgRow), ]
  # if min=max, label only 1 drug
  if (min(subjDrugRows$OgRow)==max(subjDrugRows$OgRow)){
    subjDrugRows$Drug='DrugOnly'
    subjRows$Drug[subjDrugRows$OgRow==subjDrugRows$OgRow]='Drug'
  } else {
    # Assign 'Drug' values to the original rows
    subjRows$Drug[subjRows$OgRow %in% drug1Rows$OgRow] <- 'Drug1'
    subjRows$Drug[subjRows$OgRow %in% drug2Rows$OgRow] <- 'Drug2'
  }
  ### Now pull out methyl vs. psilo
  DoseCorresp=subjDoseCorresp[subjDoseCorresp$X1==subj,]
  # if Drug1 = 1 in key, Drug1 = psilo
  if (subjDrugRows$Drug[1] == 'DrugOnly') {
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
```

    ## [1] "PS24"

``` r
# for resting-state 6
for (s in 1:length(unique(rs6m2$Subjects))){
  # this subject
  subj=unique(rs6m2$Subjects)[s]
  print(subj)
  # Subset rows for the current subject
  subjRows <- rs6m2[rs6m2$Subjects == subj, ]
  # Subset rows for the 'Drug' dosage
  subjDrugRows <- subjRows[subjRows$Dosage == 'Drug', ]
  # Identify 'Drug1' and 'Drug2' based on 'OgRow'
  drug1Rows <- subjDrugRows[subjDrugRows$OgRow == min(subjDrugRows$OgRow), ]
  drug2Rows <- subjDrugRows[subjDrugRows$OgRow == max(subjDrugRows$OgRow), ]
  # if min=max, label only 1 drug
  if (min(subjDrugRows$OgRow)==max(subjDrugRows$OgRow)){
    subjDrugRows$Drug='DrugOnly'
    subjRows$Drug[subjDrugRows$OgRow==subjDrugRows$OgRow]='Drug'
  } else {
    # Assign 'Drug' values to the original rows
    subjRows$Drug[subjRows$OgRow %in% drug1Rows$OgRow] <- 'Drug1'
    subjRows$Drug[subjRows$OgRow %in% drug2Rows$OgRow] <- 'Drug2'
  }
  ### Now pull out methyl vs. psilo
  DoseCorresp=subjDoseCorresp[subjDoseCorresp$X1==subj,]
  # if Drug1 = 1 in key, Drug1 = psilo
  if (subjDrugRows$Drug[1] == 'DrugOnly') {
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
```

    ## [1] "PS24"

``` r
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
```

``` r
### make chronological order column.
allScans$TemporalOrder=0
for (s in 1:length(unique(allScans$Subjects))){
  # get subject ID
  subject=(unique(allScans$Subjects))[s]
  allScansSubj=allScans[allScans$Subjects==subject,]
  # first logical temporal determination is in Dosage column: before < drug 1 < between < drug 2 < after
  allScansSubj$TemporalOrder[allScansSubj$Dosage=='before']=1000
  allScansSubj$TemporalOrder[allScansSubj$Dosage=='Drug']=2000
  # note drug1 vs drug 2 to be demarcated later, after this loop, using Drug column 
  allScansSubj$TemporalOrder[allScansSubj$Dosage=='between']=3000
  allScansSubj$TemporalOrder[allScansSubj$Dosage=='after']=4000
  # second logical temporal determination is OgRow column: lower og row number = earlier on
  allScansSubj$TemporalOrder <- allScansSubj$TemporalOrder + allScansSubj$OgRow
  # use both pieces of logic combined to get sequence for each subject of first to last scan (numbered 1 for first and x for last scan)
  allScansSubj$TemporalOrder <- rank(allScansSubj$TemporalOrder, ties.method = "min")
  allScansSubj_TemporalOrder<-as.numeric(allScansSubj$TemporalOrder)
  # sorted values for indexing
  allScansSubj_TemporalOrderSorted=sort(unique(allScansSubj$TemporalOrder))
  # convert to by-1 sequence
  for (i in 1:length(unique(allScansSubj$TemporalOrder))){
      numberOfInterest=allScansSubj_TemporalOrderSorted[i]
      allScansSubj_TemporalOrder[allScansSubj_TemporalOrder==numberOfInterest]=i
  }
  allScansSubj$TemporalOrder<-allScansSubj_TemporalOrder
  allScans[allScans$Subject == subject, ]$TemporalOrder <- allScansSubj$TemporalOrder
}
```

``` r
# Now we need to place drug2 scans after all betweens. Subtract # of drug2 scans from between temporal order, plop drug2 in the integer sequence space between is translated down by.


# use subjDoseCorresp X2 to find where methyl is drug1 (X2=2 indicates methyl is drug1, X2=1 indicates psilo is drug1)

# New loop to adjust TemporalOrder based on drug2 placement
for (s in 1:length(unique(allScans$Subject))) {
  # get subject ID
  subject <- (unique(allScans$Subject))[s]
  allScansSubj <- allScans[allScans$Subject == subject, ]
  # get which dose is psil/methyl
  DoseCorresp=subjDoseCorresp$X2[subjDoseCorresp$X1==subject]
  # if psilo is drug 1
  if (DoseCorresp==1){
        # recoded 6/29/24 to use same extraction as above
        # first logical temporal determination is in Dosage column: before < drug 1 < between < drug 2 < after
        allScansSubj$TemporalOrder[allScansSubj$Dosage=='before']=1000
        allScansSubj$TemporalOrder[allScansSubj$Drug=='Psilo']=2000
        # note drug1 vs drug 2 to be demarcated later, after this loop, using Drug column 
        allScansSubj$TemporalOrder[allScansSubj$Dosage=='between']=3000
        allScansSubj$TemporalOrder[allScansSubj$Drug=='Methyl']=4000
        allScansSubj$TemporalOrder[allScansSubj$Dosage=='after']=5000
        # second logical temporal determination is OgRow column: lower og row number = earlier on
        allScansSubj$TemporalOrder <- allScansSubj$TemporalOrder + allScansSubj$OgRow
        # use both pieces of logic combined to get sequence for each subject of first to last scan (numbered 1 for first and x for last scan)
        allScansSubj$TemporalOrder <- rank(allScansSubj$TemporalOrder, ties.method = "min")
        allScansSubj_TemporalOrder<-as.numeric(allScansSubj$TemporalOrder)
        # sorted values for indexing
        allScansSubj_TemporalOrderSorted=sort(unique(allScansSubj$TemporalOrder))
        # convert to by-1 sequence
        for (i in 1:length(unique(allScansSubj$TemporalOrder))){
            numberOfInterest=allScansSubj_TemporalOrderSorted[i]
            allScansSubj_TemporalOrder[allScansSubj_TemporalOrder==numberOfInterest]=i
        }
        allScansSubj$TemporalOrder<-allScansSubj_TemporalOrder
        allScans[allScans$Subject == subject, ]$TemporalOrder <- allScansSubj$TemporalOrder
     # if psilo is drug 2
  } else if (DoseCorresp==2) {
    # recoded 6/29/24 to use same extraction as above
        # first logical temporal determination is in Dosage column: before < drug 1 < between < drug 2 < after
        allScansSubj$TemporalOrder[allScansSubj$Dosage=='before']=1000
        allScansSubj$TemporalOrder[allScansSubj$Drug=='Methyl']=2000
        # note drug1 vs drug 2 to be demarcated later, after this loop, using Drug column 
        allScansSubj$TemporalOrder[allScansSubj$Dosage=='between']=3000
        allScansSubj$TemporalOrder[allScansSubj$Drug=='Psilo']=4000
        allScansSubj$TemporalOrder[allScansSubj$Dosage=='after']=5000
        # second logical temporal determination is OgRow column: lower og row number = earlier on
        allScansSubj$TemporalOrder <- allScansSubj$TemporalOrder + allScansSubj$OgRow
        # use both pieces of logic combined to get sequence for each subject of first to last scan (numbered 1 for first and x for last scan)
        allScansSubj$TemporalOrder <- rank(allScansSubj$TemporalOrder, ties.method = "min")
        allScansSubj_TemporalOrder<-as.numeric(allScansSubj$TemporalOrder)
        # sorted values for indexing
        allScansSubj_TemporalOrderSorted=sort(unique(allScansSubj$TemporalOrder))
        # convert to by-1 sequence
        for (i in 1:length(unique(allScansSubj$TemporalOrder))){
            numberOfInterest=allScansSubj_TemporalOrderSorted[i]
            allScansSubj_TemporalOrder[allScansSubj_TemporalOrder==numberOfInterest]=i
        }
        allScansSubj$TemporalOrder<-allScansSubj_TemporalOrder
        allScans[allScans$Subject == subject, ]$TemporalOrder <- allScansSubj$TemporalOrder
  }
# end loop
}
```

``` r
# need to add a session column with After1, After2, etc. to match vertexwise

# Adding session column
allScans$Session <- NA

# Unique subjects
subjects <- unique(allScans$Subjects)

# for each subject
for (subject in subjects) {
  # Subset for each subject
  subject_data <- allScans[allScans$Subjects == subject, ]
  
  # for all categories (before, between, drug, after)
  categories <- unique(subject_data$Dosage)

  # use Temporal Order column (all rows with same Temporal Order column) to set "Session" variable as concatenated category variable + 1 for first in temporal order, , cat. + 2 for second in order, etc.
  # Loop over each category
  for (category in categories) {
    # Subset for each category
    category_data <- subject_data[subject_data$Dosage == category, ]
    
    # Get unique TemporalOrder values
    temporal_orders <- unique(category_data$TemporalOrder)
    
    # Loop over each TemporalOrder
    for (i in seq_along(temporal_orders)) {
      temporal_order <- temporal_orders[i]
      # Assign session values
      category_data$Session[category_data$TemporalOrder == temporal_order] <- paste(category, i, sep = "")
    }
    
    # Update the main dataframe
    allScans[allScans$Subjects == subject & allScans$Dosage == category, ] <- category_data
  }
}

# end

# final thing is to convert all "before" to "Baseline", "between" to "Between". and "after" to "After"
allScans$Session <- gsub("before", "Baseline", allScans$Session)
allScans$Session <- gsub("between", "Between", allScans$Session)
allScans$Session <- gsub("after", "After", allScans$Session)
```

``` r
# remove data that needs to be removed (<250 TRs)
allScans=allScans[allScans$RemTRs>250,]

# make donut plot)
donutData<- data.frame(
  Category=levels(allScans$Drug)[1:3],
  count=tabulate(allScans$Drug)
)

# convert labels to be consistent across studies
donutData$Category=c('No Drug','Active Control','Psychedelic')

# percentages
donutData$frac = donutData$count / sum(donutData$count)

# Compute the cumulative for plotting
donutData$ymax = cumsum(donutData$frac)

# Compute the bottom of each rectangle (plotted as rectangle and coord_polar'ed)
donutData$ymin = c(0, head(donutData$ymax, n=-1))
 
# plot
ggplot(donutData, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Category)) +
     geom_rect() +
     coord_polar(theta="y") + 
     xlim(c(1, 4)) + theme_classic(base_size=28)+theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank()
  )+guides(fill = guide_legend(title = NULL))+scale_fill_manual(values = c("#EF9500","#002642","#840032"))
```

![](Stats_n_Viz_psil_DMNSeg_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

``` r
### match 99's to <9x ID's

# temporal order of 99's has to be after 03 16 18 19!
allScans$TemporalOrder[allScans$Subjects=='PS93']=allScans$TemporalOrder[allScans$Subjects=='PS93']+max(allScans$TemporalOrder[allScans$Subjects=='PS03'])

allScans$TemporalOrder[allScans$Subjects=='PS96']=allScans$TemporalOrder[allScans$Subjects=='PS96']+max(allScans$TemporalOrder[allScans$Subjects=='PS16'])

allScans$TemporalOrder[allScans$Subjects=='PS98']=allScans$TemporalOrder[allScans$Subjects=='PS98']+max(allScans$TemporalOrder[allScans$Subjects=='PS18'])

allScans$TemporalOrder[allScans$Subjects=='PS99']=allScans$TemporalOrder[allScans$Subjects=='PS99']+max(allScans$TemporalOrder[allScans$Subjects=='PS19'])

# this technically makes all non drug scans "after" for these PTs
allScans$Dosage[allScans$Subjects=='PS93']='after'
allScans$Dosage[allScans$Subjects=='PS96']='after'
allScans$Dosage[allScans$Subjects=='PS98']='after'
allScans$Dosage[allScans$Subjects=='PS99']='after'

# rename to formally fold into same subject ID
allScans$Subjects[allScans$Subjects=='PS93']='PS03'
allScans$Subjects[allScans$Subjects=='PS96']='PS16'
allScans$Subjects[allScans$Subjects=='PS98']='PS18'
allScans$Subjects[allScans$Subjects=='PS99']='PS19'
allScans$Subjects=as.factor(allScans$Subjects)
```

``` r
# model
fit_lme <- lme(TDProp1 ~ Drug + RemTRs + FD, random = ~ 1 | Subjects, data = allScans)
summaryLME<-summary(fit_lme)
# match to one-tailed
paste('one sided p (confirmatory of psil):', pt(summaryLME$tTable[3,4],summaryLME$tTable[3,3],lower=FALSE))
```

    ## [1] "one sided p (confirmatory of psil): 3.66529942033336e-07"

``` r
# methylphenidate vs. psilocybin
drugScans=allScans[allScans$Drug!='none',]
fit_lme <- lme(TDProp1 ~ Drug + RemTRs + FD, random = ~ 1 | Subjects, data = drugScans)
summaryLME<-summary(fit_lme)
# match to one-tailed
paste('one sided p (confirmatory of psil):', pt(summaryLME$tTable[2,4],summaryLME$tTable[2,3],lower=FALSE))
```

    ## [1] "one sided p (confirmatory of psil): 0.00318920953970199"

``` r
# save out this df for merging
saveRDS(allScans,'~/Downloads/Psil_DMNSeg_Merged.rds')
```
