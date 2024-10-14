p50 psil
================
2024-02-03

``` r
library(reshape2)
library(ggplot2)
library(visreg)
library(nlme)
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following object is masked from 'package:nlme':
    ## 
    ##     collapse

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
# prop angles
rs1=read.csv('~/Downloads/rs1_Psil_propsMerged(5).csv',header=F)
rs2=read.csv('~/Downloads/rs2_Psil_propsMerged(5).csv',header=F)
rs3=read.csv('~/Downloads/rs3_Psil_propsMerged(5).csv',header=F)
rs4=read.csv('~/Downloads/rs4_Psil_propsMerged(5).csv',header=F)
rs5=read.csv('~/Downloads/rs5_Psil_propsMerged(5).csv',header=F)
rs6=read.csv('~/Downloads/rs6_Psil_propsMerged(5).csv',header=F)
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
    ## [1] "PS21"
    ## [1] "PS24"
    ## [1] "PS93"
    ## [1] "PS96"
    ## [1] "PS98"
    ## [1] "PS99"

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
    ## [1] "PS21"
    ## [1] "PS24"
    ## [1] "PS93"
    ## [1] "PS96"
    ## [1] "PS98"
    ## [1] "PS99"
    ## [1] "PS16"

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
    ## [1] "PS99"

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
# saveout merged DF for vertexwise analyses on HPC
saveRDS(allScans,'~/forVertexwise_psil.rds')

# make donut plot)
donutData<- data.frame(
  Category=levels(allScans$Drug)[1:3],
  count=tabulate(allScans$Drug)
)

# convert labels to be consistent across studies
donutData$Category=c('Baseline','Active Control','Psychedelic')

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

![](Stats_n_Viz_psil_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

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
paste('one sided p (confirmatory of psil):', pt(summaryLME$tTable[3,4],summaryLME$tTable[3,3],lower=TRUE))
```

    ## [1] "one sided p (confirmatory of psil): 0.0162299631261041"

``` r
# methylphenidate vs. psilocybin
drugScans=allScans[allScans$Drug!='none',]
fit_lme <- lme(TDProp1 ~ Drug + RemTRs + FD, random = ~ 1 | Subjects, data = drugScans)
summaryLME<-summary(fit_lme)
# match to one-tailed
paste('one sided p (confirmatory of psil):', pt(summaryLME$tTable[2,4],summaryLME$tTable[2,3],lower=TRUE))
```

    ## [1] "one sided p (confirmatory of psil): 0.0413240200310478"

``` r
# save out this df for merging
saveRDS(allScans,'~/Downloads/Psil_BUP_Merged.rds')
```

``` r
# alternative coding to above
allScans$Condition=allScans$Dosage
# correct betweens to after where psilo has already occurred
for (s in 1:length(unique(allScans$Subject))) {
  # get subject ID
  subject <- (unique(allScans$Subject))[s]
  allScansSubj <- allScans[allScans$Subject == subject, ]
  if (min(allScansSubj$TemporalOrder[allScansSubj$Drug=='Psilo'])<min(allScansSubj$TemporalOrder[allScansSubj$Dosage=='between'])){
    allScansSubj$Dosage[allScansSubj$Dosage=='between']='after'
    # alter between visits to after if this condition is met
    allScans[allScans$Subject == subject, ]$Condition <- allScansSubj$Dosage
  }
}

# correct betweens to before where psilo has not occurred
allScans$Condition [allScans$Condition == 'between']<- 'before'

# construe after scans in terms of how many scans after
levels(allScans$Condition)<-c(levels(allScans$Condition),'immed_afterPsil','immed_afterMethyl')
# for each subj
for (s in 1:length(unique(allScans$Subject))) {
  # get subject ID
  subject <- (unique(allScans$Subject))[s]
  allScansSubj <- allScans[allScans$Subject == subject, ]
  TempOrdersubj=allScansSubj$TemporalOrder
  # find psil scan
  psilScanInds=TempOrdersubj[allScansSubj$Drug=='Psilo']
  unqPsilScans=unique(psilScanInds)
  for (unqPsilScans in 1:length(unqPsilScans)){
      # psil scan temporal index
      p_s_t_i=psilScanInds[unqPsilScans]
      # mark after scan as immedafter in cond
      allScansSubj$Condition[allScansSubj$TemporalOrder==(p_s_t_i+1)]='immed_afterPsil'
      allScansSubj$Condition[allScansSubj$TemporalOrder==(p_s_t_i+2)]='immed_afterPsil'
  }
  # find methyl scan
  methScanInds=TempOrdersubj[allScansSubj$Drug=='Methyl']
  unqmethScans=unique(methScanInds)
  for (unqmethScans in 1:length(unqmethScans)){
      # meth scan temporal index
      m_s_t_i=methScanInds[unqmethScans]
      # mark after scan as immedafter in cond
      allScansSubj$Condition[allScansSubj$TemporalOrder==(m_s_t_i+1)]='immed_afterMethyl'
      allScansSubj$Condition[allScansSubj$TemporalOrder==(m_s_t_i+2)]='immed_afterMethyl'
  }
  # insert back into master dataframe
  allScans[allScans$Subject == subject, ]$Condition <- allScansSubj$Condition
}
  

# after condition no longer meaningfully distinct from before
allScans$Condition2<-as.factor(allScans$Condition)
allScans$Condition2[allScans$Condition2=='after']='before'
levels(allScans$Condition2)=c(levels(allScans$Condition2),'Non-Drug','Psilo','Methyl')
allScans$Condition2[allScans$Condition2=='before']='Non-Drug'
allScans$Condition2[allScans$Drug=='Methyl']='Methyl'
allScans$Condition2[allScans$Drug=='Psilo']='Psilo'
```

``` r
# model it with non-drug as reference scans
allScans <- within(allScans, Condition2 <- relevel(Condition2, ref = 6))
```

``` r
# figure 2 plots: drug vs. nondrug
allScans$Psilo=0
allScans$Psilo[allScans$Drug=='Psilo']=1
allScans$Psilo<-as.factor(allScans$Psilo)


# generate extended color pal for subject plotting
library(grDevices)
# Define the extended custom palette function
extended_palette <- colorRampPalette(rev(c("#FFEE00", "#EF9500", "#002642", "#c1004f", "#000000")))

# Generate a palette with the number of unique levels in V1
num_colors <- length(unique(allScans$Subjects))
generated_colors <- extended_palette(num_colors)

# get unique subj names
allScans$Subjects <- droplevels(allScans$Subjects)
unique_values <- unique(allScans$Subjects)
new_labels <- paste0("PSIL", seq_along(unique_values))
names(new_labels) <- unique_values

# Replace the values in Subjects using the new labels
allScans$People <- new_labels[allScans$Subjects]

# final ordering
allScans$People <- factor(allScans$People, levels = new_labels)

# just non-drug for easier interpretability for this plot
allScansNoMeth=allScans[allScans$Drug!='Methyl',]
# plot residuals from reduced model: lme(TDProp1 ~ Drug + RemTRs + FD)
reducedModel=lme(TDProp1 ~RemTRs + FD, random = ~ 1 | Subjects, data = allScansNoMeth)
allScansNoMeth$ResidualsReduced=resid(reducedModel)+mean(allScansNoMeth$TDProp1)

# gaussian jitter
set.seed(1)
allScansNoMeth$JitteredPsilo <- as.numeric(allScansNoMeth$Psilo) + rnorm(nrow(allScansNoMeth), mean = 0, sd = 0.1)

ggplot(allScansNoMeth, aes(x = JitteredPsilo, y = ResidualsReduced)) +
  geom_point(alpha = 0.8, size = 4, aes(color = People)) +  # Points with Gaussian jitter
  geom_boxplot(aes(group = Psilo), alpha = 0.2, outlier.shape = NA, width = 0.25) +  # Boxplot
  labs(title = "Psilocybin vs. No-drug scans \n",
       x = "",
       y = "% Bottom-up") + 
  scale_x_continuous(breaks = 1:2, labels = c('No Drug','Psilocybin')) +
  theme_minimal(base_size=26)+scale_color_manual(values = generated_colors)
```

![](Stats_n_Viz_psil_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

``` r
# job application version - 350 x 600
ggplot(allScansNoMeth, aes(x = Psilo, y = ResidualsReduced)) +
  geom_jitter(width = 0.25, height = 0, alpha = 0.8, size = 2, aes(color = People)) +  # Jittered points
  geom_boxplot(alpha = 0.2, outlier.shape = NA) +     # Boxplot
  labs(x = "",
       y = "% Bottom-up") +
  scale_x_discrete(labels = c('No Drug', 'Psilocybin')) +
  scale_color_manual(values = generated_colors) +  # Custom generated color palette
  theme_minimal(base_size = 28)+
  theme(legend.position = "none",axis.text.x=element_text(angle=45))
```

![](Stats_n_Viz_psil_files/figure-gfm/unnamed-chunk-15-2.png)<!-- -->

``` r
# full model for stats
fit_lme <- lme(TDProp1 ~ Psilo + RemTRs + FD, random = ~ 1 | Subjects, data = allScansNoMeth)
summaryLME<-summary(fit_lme)
# match to one-tailed
paste('one sided p (confirmatory of psil):', pt(summaryLME$tTable[2,4],summaryLME$tTable[2,3],lower=TRUE))
```

    ## [1] "one sided p (confirmatory of psil): 0.0202306673445177"

``` r
# figure 2 plots: methylphenidate vs. psilocybin

# subset out sessions where participant received pill
drugScans=allScans[allScans$Drug!='none',]

# plot residuals from reduced model: lme(TDProp1 ~ Drug + RemTRs + FD)
reducedModel=lme(TDProp1 ~  RemTRs + FD, random = ~ 1 | Subjects, data = drugScans)
drugScans$ResidualsReduced=resid(reducedModel)+mean(drugScans$TDProp1)

# gaussian jitter, as.numeric - 1 because none is retained as factor level despite not being populated
set.seed(1)
drugScans$JitteredDrug <- as.numeric(drugScans$Drug)-1 + rnorm(nrow(drugScans), mean = 0, sd = 0.1)

ggplot(drugScans, aes(x = JitteredDrug, y = ResidualsReduced)) +
  geom_point(alpha = 0.8, size = 4, aes(color = People)) +  # Points with Gaussian jitter
  geom_boxplot(aes(group = Drug), alpha = 0.2, outlier.shape = NA, width = 0.25) +  # Boxplot
  labs(title = "Psilocybin vs. Active Placebo \n \n",
       x = "",
       y = "% Bottom-up") + 
  scale_x_continuous(breaks = 1:2, labels = c('Methylphenidate','Psilocybin')) +
  theme_minimal(base_size=26)+scale_color_manual(values = generated_colors)+theme(axis.text.x = element_text(size = 20))
```

![](Stats_n_Viz_psil_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

``` r
# full model for stats
fit_lme <- lme(TDProp1 ~ Drug + RemTRs + FD, random = ~ 1 | Subjects, data = drugScans)
summaryLME<-summary(fit_lme)
# match to one-tailed
paste('one sided p (confirmatory of psil):', pt(summaryLME$tTable[2,4],summaryLME$tTable[2,3],lower=TRUE))
```

    ## [1] "one sided p (confirmatory of psil): 0.0413240200310478"

``` r
# figure 2 plots: follow-up scans
# get residuals
reducedModel=lme(TDProp1 ~RemTRs + FD, random = ~ 1 | Subjects, data = allScans)
allScans$ResidualsReduced=resid(reducedModel)+mean(allScans$TDProp1)
# plot one version without methyl conditions for main text
allScansNoMeth=allScans[allScans$Condition2!='Methyl',]
allScansNoMeth=allScansNoMeth[allScansNoMeth$Condition2!='immed_afterMethyl',]

ggplot(allScansNoMeth, aes(x = Condition2, y = ResidualsReduced)) +
  geom_jitter(width = 0.25, height = 0, alpha = 0.8,size=3,aes(color=People)) +  # Jittered points
  geom_boxplot(alpha = 0.2,outlier.shape = NA) +     # Boxplot
  labs(title = "Psilocybin: Lasting effects \n",
       x = "",
       y = "% Bottom-up") + scale_x_discrete(labels=c('No Drug','Post-Psil','During Psil'))+
  theme_minimal(base_size=25)+scale_color_manual(values = generated_colors)
```

![](Stats_n_Viz_psil_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

``` r
# full model for stats
fit_lme <- lme(TDProp1 ~ Condition2+RemTRs+FD, random = ~ 1 | Subjects, data = allScans)
summaryLME<-summary(fit_lme)
# match to one-tailed
paste('one sided p (confirmatory of psil):', pt(summaryLME$tTable[4,4],summaryLME$tTable[4,3],lower=TRUE))
```

    ## [1] "one sided p (confirmatory of psil): 0.00318557373815308"

``` r
paste('one sided p (confirmatory of post-psil):', pt(summaryLME$tTable[3,4],summaryLME$tTable[3,3],lower=TRUE))
```

    ## [1] "one sided p (confirmatory of post-psil): 0.0202128684948413"

``` r
# figure 3: load in all DMN measures
BUP=readRDS('~/Downloads/Psil_BUP_Merged.rds')
DMNSeg=readRDS('~/Downloads/Psil_DMNSeg_Merged.rds')
Mag=readRDS('~/Downloads/Psil_MagMerged.rds')
TA=readRDS('~/Downloads/Psil_TAuto_Merged.rds')
# merge em all
colnames(BUP)[1]<-'BUP'
colnames(DMNSeg)[1]<-'DMNFC'
colnames(Mag)[1]<-'Magnitudes'
colnames(TA)[1]<-'TmpAutoCor'

# merge em
PsilMerged=merge(BUP,DMNSeg,by=c('RemTRs','OgRow','FD','Subjects','Task','Dosage','TemporalOrder','Session','Drug'))
PsilMerged=merge(PsilMerged,Mag,by=c('RemTRs','OgRow','FD','Subjects','Task','Dosage','TemporalOrder','Session','Drug'))
PsilMerged=merge(PsilMerged,TA,by=c('RemTRs','OgRow','FD','Subjects','Task','Dosage','TemporalOrder','Session','Drug'))
```

    ## Warning in merge.data.frame(PsilMerged, TA, by = c("RemTRs", "OgRow", "FD", :
    ## column names 'TDProp2.x', 'TDProp3.x', 'TDProp4.x', 'TDProp2.y', 'TDProp3.y',
    ## 'TDProp4.y' are duplicated in the result

``` r
### Psil vs. methyl version
PsilMergedDrug=PsilMerged[PsilMerged$Drug!='none',]

# inclusive models
fit_BUP <- lme(BUP ~ Drug + RemTRs + FD + Magnitudes, random = ~ 1 | Subjects, data = PsilMergedDrug)
fit_Mag <- lme(Magnitudes ~ Drug + RemTRs + FD + BUP, random = ~ 1 | Subjects, data = PsilMergedDrug)

### Psil vs. no drug version
# remove methylphenidate for ease of interpretation
PsilMerged=PsilMerged[PsilMerged$Drug!='Methyl',]

fit_BUP <- lme(BUP ~ Drug + RemTRs + FD + Magnitudes, random = ~ 1 | Subjects, data = PsilMerged)
fit_Mag <- lme(Magnitudes ~ Drug + RemTRs + FD + BUP, random = ~ 1 | Subjects, data = PsilMerged)


PsilMerged$DrugBinary=0
PsilMerged$DrugBinary[PsilMerged$Drug=='Psilo']=1



# AUC thangs
library(pROC)
```

    ## Type 'citation("pROC")' for a citation.

    ## 
    ## Attaching package: 'pROC'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     cov, smooth, var

``` r
library(plotROC)
```

    ## 
    ## Attaching package: 'plotROC'

    ## The following object is masked from 'package:pROC':
    ## 
    ##     ggroc

``` r
# Fit logistic regression models
model1 <- glm(DrugBinary ~ FD + RemTRs + DMNFC+TmpAutoCor, data = PsilMerged, family = binomial)
model2 <- glm(DrugBinary ~ FD + RemTRs + DMNFC + BUP + TmpAutoCor+Magnitudes, data = PsilMerged, family = binomial)

# Predict probabilities
prob1 <- predict(model1, type = "response")
prob2 <- predict(model2, type = "response")
```

``` r
# Calculate AUC for each model
roc1 <- roc(PsilMerged$DrugBinary, prob1)
```

    ## Setting levels: control = 0, case = 1

    ## Setting direction: controls < cases

``` r
roc2 <- roc(PsilMerged$DrugBinary, prob2)
```

    ## Setting levels: control = 0, case = 1
    ## Setting direction: controls < cases

``` r
# Print AUC values
auc1 <- auc(roc1)
auc2 <- auc(roc2)

print(paste("AUC for DMN Correlations:", auc1))
```

    ## [1] "AUC for DMN Correlations: 0.767123287671234"

``` r
print(paste("AUC for DMN Propagations:", auc2))
```

    ## [1] "AUC for DMN Propagations: 0.888698630136984"

``` r
# Calculate AUC difference between full and reduced models
auc_diff <- auc2 - auc1

# Create a combined data frame for all models
df <- data.frame(
  labels = as.numeric(rep(PsilMerged$Drug, 2)),
  predictions = c(prob1, prob2),
  model = factor(rep(c("DMN Correlations", "+DMN Propagations"), each = nrow(PsilMerged)))
)

# Generate the ROC plot
ggplot(df, aes(m = predictions, d = labels, color = model)) + 
  geom_roc(n.cuts = 0, labels = FALSE) + 
  ylim(0, 1) + ylab('True Positive Rate') +xlab('False Positive Rate')+
  ggtitle("ROC Curves for Predicting Psilocybin") + 
  theme_minimal(base_size=18) + 
  scale_color_manual(values = c("#09416b","#c12139"))+
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray")+
  theme(legend.position = "none")
```

    ## Warning in verify_d(data$d): D not labeled 0/1, assuming 1 = 0 and 3 = 1!

![](Stats_n_Viz_psil_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->

``` r
# make equivalent AUC calculations on permuted data

#### commented out because it takes forecer to run, but does work when uncommented

## initialize AUC difference vectors
#auc_diffs <- rep(NA, 1000)
#
## 1. permute each DMN variable (and FD)
#set.seed(1)
## temp to allow this to run quick, analysis ran on 1000
#for (i in 1:2){
##for (i in 1:1000){
#  print(i)
#  # permute DMNMag
#  PsilMerged$DMNMag_perm <- sample(PsilMerged$Magnitudes)
#  # permute TDProp1
#  PsilMerged$TDProp1_perm <- sample(PsilMerged$BUP)
# 
#  # Fit logistic regression models
#  model1 <- glm(DrugBinary ~ FD + RemTRs + DMNFC+TmpAutoCor, data = PsilMerged, family = binomial)
#  model2_perm <- glm(DrugBinary ~ FD + RemTRs + DMNFC+TmpAutoCor+TDProp1_perm+DMNMag_perm, data = PsilMerged, family = binomial)
#   
#  # 3. calculate AUC difference between full and reduced models with permuted data
#  roc1 <- roc(PsilMerged$DrugBinary, predict(model1, type = "response"))
#  roc2_perm <- roc(PsilMerged$DrugBinary, predict(model2_perm, type = "response"))
#  
#  # Print AUC values
#  auc1 <- auc(roc1)
#  auc2_perm <-auc(roc2_perm)
#
#  # populate auc_diff vectors
#  # DMN correlations vs. full (permuted props) model
#  auc_diffs[i] <- auc1 - auc2_perm
#}
## 4. Compare true AUC differences to permuted AUC differences
#
#sum(auc_diffs>auc_diff)
# 0 indicates p <0.001
```

``` r
# full data models
ta_model <- lme(TmpAutoCor ~ FD + Drug+RemTRs, random = ~ 1 | Subjects, data = PsilMerged)
ds_model <- lme(DMNFC ~ FD + Drug+RemTRs, random = ~ 1 | Subjects, data = PsilMerged)
```

``` r
# and bootstrapping: also figure 3
# Set the number of bootstrap samples
num_bootstrap_samples <- 1000
# initialize t vectors
td_d<-rep(0,num_bootstrap_samples)
td_fd<-rep(0,num_bootstrap_samples)
ta_d<-rep(0,num_bootstrap_samples)
ta_fd<-rep(0,num_bootstrap_samples)
cx_d<-rep(0,num_bootstrap_samples)
cx_fd<-rep(0,num_bootstrap_samples)
ds_d<-rep(0,num_bootstrap_samples)
ds_fd<-rep(0,num_bootstrap_samples)
dm_d<-rep(0,num_bootstrap_samples)
dm_fd<-rep(0,num_bootstrap_samples)

# bootstrap loops
set.seed(1)
for (i in 1:num_bootstrap_samples){
  # resample data
  data=PsilMerged[sample(nrow(PsilMerged), replace = TRUE), ]
  # fit on all models
  td_model <- lme(BUP ~ FD + Drug+RemTRs, random = ~ 1 | Subjects, data = data)
  ta_model <- lme(TmpAutoCor ~ FD + Drug+RemTRs, random = ~ 1 | Subjects, data = data)
  ds_model <- lme(DMNFC ~ FD + Drug+RemTRs, random = ~ 1 | Subjects, data = data)
  dm_model <- lme(Magnitudes ~ FD + Drug+RemTRs, random = ~ 1 | Subjects, data = data)
  # get t-values
  td_d[i]=summary(td_model)$tTable[rownames(summary(td_model)$tTable) == "DrugPsilo", "t-value"]
  td_fd[i]=summary(td_model)$tTable[rownames(summary(td_model)$tTable) == "FD", "t-value"]
  ta_d[i]=summary(ta_model)$tTable[rownames(summary(ta_model)$tTable) == "DrugPsilo", "t-value"]
  ta_fd[i]=summary(ta_model)$tTable[rownames(summary(ta_model)$tTable) == "FD", "t-value"]
  ds_d[i]=summary(ds_model)$tTable[rownames(summary(ds_model)$tTable) == "DrugPsilo", "t-value"]
  ds_fd[i]=summary(ds_model)$tTable[rownames(summary(ds_model)$tTable) == "FD", "t-value"]
  dm_d[i]=summary(dm_model)$tTable[rownames(summary(dm_model)$tTable) == "DrugPsilo", "t-value"]
  dm_fd[i]=summary(dm_model)$tTable[rownames(summary(dm_model)$tTable) == "FD", "t-value"]
}
# convert to dataframes
td_d=data.frame(td_d)
ta_d=data.frame(ta_d)
ds_d=data.frame(ds_d)
dm_d=data.frame(dm_d)
td_fd=data.frame(td_fd)
ta_fd=data.frame(ta_fd)
ds_fd=data.frame(ds_fd)
dm_fd=data.frame(dm_fd)

colnames(td_d)='tstat'
colnames(ta_d)='tstat'
colnames(ds_d)='tstat'
colnames(dm_d)='tstat'
colnames(td_fd)='tstat'
colnames(ta_fd)='tstat'
colnames(ds_fd)='tstat'
colnames(dm_fd)='tstat'

# set column names for merging
td_d$Cov='Drug'
ta_d$Cov='Drug'
ds_d$Cov='Drug'
dm_d$Cov='Drug'

td_fd$Cov='FD'
ta_fd$Cov='FD'
ds_fd$Cov='FD'
dm_fd$Cov='FD'

td_d$Model='Bottom-up %'
ta_d$Model='AutoCor'
ds_d$Model='Integration'
dm_d$Model='Magnitude'

td_fd$Model='Bottom-up %'
ta_fd$Model='AutoCor'
ds_fd$Model='Integration'
dm_fd$Model='Magnitude'

bootstrap_results_FD=rbind(td_fd,ta_fd,ds_fd,dm_fd)
bootstrap_results_Drug=rbind(td_d,ta_d,ds_d,dm_d)
# Calculate the average t-value for each Model category
average_t_values <- bootstrap_results_FD %>%
  group_by(Model) %>%
  summarize(avg_t_value = mean(tstat, na.rm = TRUE))

# Reorder the Model factor based on the average t-values
bootstrap_results_FD$Model <- factor(bootstrap_results_FD$Model, 
                                     levels = c('Integration','Bottom-up %','AutoCor','Magnitude'))

# Generate the plot
library(ggdist)
ggplot(bootstrap_results_FD, aes(x = Model, y = tstat, fill = Cov)) +
    geom_boxplot(
        width = 0.12,
        # Removing outliers
        outlier.color = NA,
        fill='#EF9500'
    ) +
    stat_dots(
        # Plotting on left side
        side = "left",
        # Adjusting position
        justification = 1.1,
        # Adjust grouping (binning) of observations
        binwidth = 0.08
    ) +
    labs(x = "Model", y = "T-Values", title = "Bootstrap T-Values for FD effect") +
    theme_minimal(base_size = 18) +
    geom_hline(yintercept = 0, linetype = "dashed")+
    # just to prevent extra x-axis expansion
    coord_cartesian(xlim = c(.8, length(unique(bootstrap_results_Drug$Model))))+
    theme(legend.position = "none")
```

    ## Warning: The provided binwidth will cause dots to overflow the boundaries of the
    ## geometry.
    ## → Set `binwidth = NA` to automatically determine a binwidth that ensures dots
    ##   fit within the bounds,
    ## → OR set `overflow = "compress"` to automatically reduce the spacing between
    ##   dots to ensure the dots fit within the bounds,
    ## → OR set `overflow = "keep"` to allow dots to overflow the bounds of the
    ##   geometry without producing a warning.
    ## ℹ For more information, see the documentation of the `binwidth` and `overflow`
    ##   arguments of `?ggdist::geom_dots()` or the section on constraining dot sizes
    ##   in vignette("dotsinterval") (`vignette(ggdist::dotsinterval)`).

![](Stats_n_Viz_psil_files/figure-gfm/unnamed-chunk-22-1.png)<!-- -->

``` r
# Calculate the average t-value for each Model category
average_t_values <- bootstrap_results_Drug %>%
  group_by(Model) %>%
  summarize(avg_t_value = mean(tstat, na.rm = TRUE))

# Reorder the Model factor based on the average t-values
bootstrap_results_Drug$Model <- factor(bootstrap_results_Drug$Model, 
                                     levels = c('Integration','AutoCor','Bottom-up %','Magnitude'))

# add fill column
bootstrap_results_Drug$Fill='DMN Correlatons'
bootstrap_results_Drug$Fill[bootstrap_results_Drug$Model=='Bottom-up %']='DMN Propagations'
bootstrap_results_Drug$Fill[bootstrap_results_Drug$Model=='Magnitude']='DMN Propagations'


# Generate the plot
ggplot(bootstrap_results_Drug, aes(x = Model, y = tstat, fill = Fill)) +
    geom_boxplot(
        width = 0.12,
        # Removing outliers
        outlier.color = NA) +
        scale_fill_manual(values=c("#c12139","#09416b"))+
    stat_dots(
        # Plotting on left side
        side = "left",
        # Adjusting position
        justification = 1.1,
        # Adjust grouping (binning) of observations
        binwidth = 0.08,
        overflow = "compress"
    ) +
    labs(x = "Model", y = "T-Values", title = "Bootstrap T-Values for Psilocybin effect") +
    theme_minimal(base_size = 18) +
    geom_hline(yintercept = 0, linetype = "dashed")+
    # just to prevent extra x-axis expansion
    coord_cartesian(xlim = c(1, length(unique(bootstrap_results_Drug$Model))))+
    theme(legend.position = "none")+ylim(c(-10.5,10.5))
```

![](Stats_n_Viz_psil_files/figure-gfm/unnamed-chunk-22-2.png)<!-- -->