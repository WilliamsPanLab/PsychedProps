MDMA main
================
2023-07-13

``` r
library(reshape2)
library(ggplot2)
library(visreg)
library(nlme)
```

``` r
# prop angles
rs1=read.csv('~/Downloads/rs1_propsMerged.csv',header=F)
rs2=read.csv('~/Downloads/rs2_propsMerged.csv',header=F)
emo=read.csv('~/Downloads/emotion_propsMerged.csv',header=F)
gambling=read.csv('~/Downloads/gambling_propsMerged.csv',header=F)
wm=read.csv('~/Downloads/wm_propsMerged.csv',header=F)

rs1$Task='rs'
rs2$Task='rs2'
emo$Task='emotion'
gambling$Task='gambling'
wm$Task='wm'

# remove subj 4
rs1=rs1[-c(4),]
rs2=rs2[-c(4),]
emo=emo[-c(4),]
gambling=gambling[-c(4),]
wm=wm[-c(4),]

# this is going to be ugly but simple
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

emobv=data.frame(cbind(emo$V1,emo$V2,emo$V3,emo$V4,emo$V17))
emop=data.frame(cbind(emo$V5,emo$V6,emo$V7,emo$V8,emo$V18))
emom1=data.frame(cbind(emo$V9,emo$V10,emo$V11,emo$V12,emo$V19))
emom2=data.frame(cbind(emo$V13,emo$V14,emo$V15,emo$V16,emo$V20))
colnames(emobv)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs')
colnames(emop)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs')
colnames(emom1)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs')
colnames(emom2)=c('TDProp1','TDProp2','TDProp3','TDProp4','RemTRs')

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

emobv$Subjects=alff$SubjID
emop$Subjects=alff$SubjID
emom1$Subjects=alff$SubjID
emom2$Subjects=alff$SubjID

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

emobv$Task='emotion'
emop$Task='emotion'
emom1$Task='emotion'
emom2$Task='emotion'

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

emobv$Dosage='baseline'
emop$Dosage='Placebo'
emom1$Dosage='80mg'
emom2$Dosage='120mg'

gamblingbv$Dosage='baseline'
gamblingp$Dosage='Placebo'
gamblingm1$Dosage='80mg'
gamblingm2$Dosage='120mg'

wmbv$Dosage='baseline'
wmp$Dosage='Placebo'
wmm1$Dosage='80mg'
wmm2$Dosage='120mg'

# combine all
allScans=rbind(rs1bv,rs1p,rs1m1,rs1m2,rs2bv,rs2p,rs2m1,rs2m2,emobv,emop,emom1,emom2,gamblingbv,gamblingp,gamblingm1,gamblingm2,wmbv,wmp,wmm1,wmm2)

# read in motion
mot=read.csv('~/Desktop/MDMA_spikes_summary.csv')

# motion merge
mergedDf=merge(mot,allScans,by=c('Subjects','Task','Dosage'))
mergedDf$Drug=0
mergedDf$Drug[mergedDf$Dosage=="120mg"]=1
mergedDf$Drug[mergedDf$Dosage=="80mg"]=1
mergedDf$Drug=as.factor(mergedDf$Drug)
mergedDf$Subjects<-as.factor(mergedDf$Subjects)
mergedDf$Dosage<-as.factor(mergedDf$Dosage)
mergedDfProps=mergedDf
```

``` r
# remove data that needs to be removed (subs 6 and 10, <250 TRs)
mergedDf=mergedDf[mergedDf$Subjects!='sub-MDMA006',]
mergedDf=mergedDf[mergedDf$Subjects!='sub-MDMA010',]
mergedDf=mergedDf[mergedDf$RemTRs>250,]
```

``` r
# assign clearer subject labels
# get unique subj names
mergedDf$Subjects <- droplevels(mergedDf$Subjects)
unique_values <- unique(mergedDf$Subjects)
new_labels <- paste0("MDMA", seq_along(unique_values))
names(new_labels) <- unique_values

# Replace the values in Subjects using the new labels
mergedDf$People <- new_labels[mergedDf$Subjects]

# final ordering
mergedDf$People <- factor(mergedDf$People, levels = new_labels)

# change rs2 to rs for accurate task-modeling
mergedDf$Task[mergedDf$Task=='rs2']='rs'
mergedDf$Task=as.factor(mergedDf$Task)
mergedDf_clean=mergedDf[mergedDf$Dosage!='baseline',]

# make donut plot
donutData<- data.frame(
  Category=levels(mergedDf$Dosage),
  count=tabulate(mergedDf$Dosage)
)

# convert labels to be consistent across studies
donutData$Category=c('Psychedelic','x','Baseline','Placebo')
# merge 120 and 80mg for clear plots
donutData$count[donutData$Category=='Psychedelic']=donutData$count[donutData$Category=='Psychedelic']+donutData$count[donutData$Category=='x']
donutData=donutData[donutData$Category!='x',]

# percentages
donutData$frac = donutData$count / sum(donutData$count)

# Compute the cumulative for plotting
donutData$ymax = cumsum(donutData$frac)

# Compute the bottom of each rectangle (plotted as rectangle and coord_polar'ed)
donutData$ymin = c(0, head(donutData$ymax, n=-1))

# convert to factor for manual ordering
donutData$Category <- factor(donutData$Category, levels = c('Placebo', 'Baseline', 'Psychedelic'))


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

![](Stats_n_Viz_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
# generate extended color pal for subject plotting
library(grDevices)
# Define the extended custom palette function
extended_palette <- colorRampPalette(rev(c("#FFEE00", "#EF9500", "#002642", "#c1004f", "#000000")))

# Generate a palette with the number of unique levels in V1
num_colors <- length(unique(mergedDf_clean$Subjects))
generated_colors <- extended_palette(num_colors)

# get head-motion regressed values for plotting
model_to_reg <- lm(TDProp1 ~ MeanFD + RemTRs+Task, data = mergedDf_clean)
mergedDf_clean$Residuals<-resid(model_to_reg)+mean(mergedDf_clean$TDProp1)

# updated with manual Gaussian jitter
# Add Gaussian jitter to x values (e.g., Drug factor levels)
set.seed(2)
mergedDf_clean$JitteredDrug <- as.numeric(mergedDf_clean$Drug) + rnorm(nrow(mergedDf_clean), mean = 0, sd = 0.1)

# figure 2a
ggplot(mergedDf_clean, aes(x = JitteredDrug, y = Residuals)) +
  geom_point(alpha = 0.8, size = 4, aes(color = People)) +  # Points with Gaussian jitter
  geom_boxplot(aes(group = Drug), alpha = 0.2, outlier.shape = NA, width = 0.25) +  # Boxplot
  labs(title = "MDMA vs. Placebo \n",
       x = "",
       y = "% Bottom-up") + 
  scale_x_continuous(breaks = 1:2, labels = c('Placebo', 'MDMA')) +
  scale_color_manual(values = generated_colors)+
  theme_minimal(base_size=28)
```

![](Stats_n_Viz_files/figure-gfm/unnamed-chunk-4-2.png)<!-- -->

``` r
# job application version - 350 x 650
ggplot(mergedDf_clean, aes(x = Drug, y = Residuals)) +
  geom_jitter(width = 0.25, height = 0, alpha = 0.8, size = 2, aes(color = People)) +  # Jittered points
  geom_boxplot(alpha = 0.2, outlier.shape = NA) +     # Boxplot
  labs(x = "",
       y = "% Bottom-up") +
  scale_x_discrete(labels = c('Placebo', 'MDMA')) +
  scale_color_manual(values = generated_colors) +  # Custom generated color palette
  theme_minimal(base_size = 28)+
  theme(legend.position = "none",axis.text.x=element_text(angle=45))
```

![](Stats_n_Viz_files/figure-gfm/unnamed-chunk-4-3.png)<!-- -->

``` r
# full model for stats
fit_lme <- lme(TDProp1 ~ Drug + RemTRs + MeanFD+Task, random = ~ 1 | Subjects, data = mergedDf_clean)
summaryLME<-summary(fit_lme)

# testing lme4 for robustness
library(lme4)
```

    ## Loading required package: Matrix

    ## 
    ## Attaching package: 'lme4'

    ## The following object is masked from 'package:nlme':
    ## 
    ##     lmList

``` r
library(lmerTest)
```

    ## 
    ## Attaching package: 'lmerTest'

    ## The following object is masked from 'package:lme4':
    ## 
    ##     lmer

    ## The following object is masked from 'package:stats':
    ## 
    ##     step

``` r
fit_lmer <- lmer(TDProp1 ~ Drug + RemTRs + MeanFD+Task + (1 | Subjects), data = mergedDf_clean)
```

    ## Warning: Some predictor variables are on very different scales: consider
    ## rescaling
    ## Warning: Some predictor variables are on very different scales: consider
    ## rescaling

``` r
# checks out

#####
#### with all non-drug scans: Figure 2b
#####

# get head-motion regressed values for plotting
model_to_reg <- lm(TDProp1 ~ MeanFD + RemTRs+Task, data = mergedDf)
mergedDf$Residuals<-resid(model_to_reg)+mean(mergedDf$TDProp1)

# updated with manual Gaussian jitter
# Add Gaussian jitter to x values (e.g., Drug factor levels)
set.seed(1)
mergedDf$JitteredDrug <- as.numeric(mergedDf$Drug) + rnorm(nrow(mergedDf), mean = 0, sd = 0.1)

# figure 2a
ggplot(mergedDf, aes(x = JitteredDrug, y = Residuals)) +
  geom_point(alpha = 0.8, size = 4, aes(color = People)) +  # Points with Gaussian jitter
  geom_boxplot(aes(group = Drug), alpha = 0.2, outlier.shape = NA, width = 0.25) +  # Boxplot
  labs(title = "MDMA vs. No-drug scans \n",
       x = "",
       y = "% Bottom-up") + 
  scale_x_continuous(breaks = 1:2, labels = c('No Drug', 'MDMA')) +
  theme_minimal(base_size = 25) +
  scale_color_manual(values = generated_colors)
```

![](Stats_n_Viz_files/figure-gfm/unnamed-chunk-4-4.png)<!-- -->

``` r
# figure 2a
#ggplot(mergedDf, aes(x = Drug, y = Residuals)) +
#  geom_jitter(width = 0.25, height = 0, alpha = 0.8,size=4,aes(color = People)) +  # Jittered points
#  geom_boxplot(alpha = 0.2,outlier.shape = NA) +     # Boxplot
#  labs(title = "MDMA vs. No-drug scans \n",
#       x = "",
#       y = "% Bottom-up") + scale_x_discrete(labels=c('No Drug','MDMA'))+
#  theme_minimal(base_size=25)+scale_color_manual(values = generated_colors)

# full model for stats
fit_lme <- lme(TDProp1 ~ Drug + RemTRs + MeanFD+Task, random = ~ 1 | Subjects, data = mergedDf)
summaryLME<-summary(fit_lme)

# lme4 test for robustness
fit_lmer <- lmer(TDProp1 ~ Drug + RemTRs + MeanFD+Task + (1 | Subjects), data = mergedDf)
```

    ## Warning: Some predictor variables are on very different scales: consider
    ## rescaling
    ## Warning: Some predictor variables are on very different scales: consider
    ## rescaling

``` r
# checks out

# make a "people"-"subjects" equivalence data frame to reference in figures down below
unique_pairs <- unique(mergedDf[c("People", "Subjects")])
```

``` r
# extract standout sessions/participants
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
# Calculate the difference in TDProp1 between Drug and Placebo conditions
df_diff <- mergedDf_clean %>%
  group_by(Subjects) %>%
  summarise(Diff_TDProp1 = mean(TDProp1[Drug == 0]-mean(TDProp1[Drug == 1]), na.rm = TRUE))

# Identify the participant with the greatest increase in TDProp1
max_diff_subject <- df_diff %>%
  filter(Diff_TDProp1 == max(Diff_TDProp1, na.rm = TRUE)) %>%
  pull(Subjects)

print(max_diff_subject)
```

    ## [1] sub-MDMA017
    ## 14 Levels: sub-MDMA001 sub-MDMA002 sub-MDMA003 sub-MDMA005 ... sub-MDMA017

``` r
# added chunk to code which sessions are post-mdma in macro timeline of study

# initialize new column in mergedDf
mergedDf$PostMDMASession <- 2

# load in subjSeshDoseCorrep (read delim, sep=' ')
subjSeshDoseCorrep <- read.table("~/Downloads/subjSeshDoseCorresp.csv", sep = ' ')

# for each subject
for (i in 1:nrow(subjSeshDoseCorrep)) {
  # if i <10
  if(i<10){
    subjectID <- paste0("sub-MDMA00", i)  # Reconstruct subject ID
  } else if(i>9){
    subjectID <- paste0("sub-MDMA0", i)  # Reconstruct subject ID
  }
  # get position of placebo
  placeboPos <- which(subjSeshDoseCorrep[i, 2:5] == "ses-01")
  print(placeboPos)
  # get position in mergedDf
  placeboRows <- mergedDf$Subject == subjectID & mergedDf$Dosage == "Placebo"
  # baseline is pos one, placebo is pos 2, 80 is pos 3, 120 is pos 4
  # if placebo is ses-01, leave as initialize 0
  # if placebo is ses-02, 1 for post-drugMacro
  # if placebo is ses-03, 1 for post-drugMacro
  if (placeboPos == 2) {
      mergedDf$PostMDMASession[placeboRows] <- 0
    } else if (placeboPos == 3 ) {
      mergedDf$PostMDMASession[placeboRows] <- 1
    } else if (placeboPos == 4) {
      mergedDf$PostMDMASession[placeboRows] <- 1
    }
}
```

    ## [1] 3
    ## [1] 2
    ## [1] 3
    ## [1] 2
    ## [1] 4
    ## [1] 3
    ## [1] 3
    ## [1] 2
    ## [1] 2
    ## [1] 4
    ## [1] 2
    ## [1] 3
    ## [1] 4
    ## [1] 4
    ## [1] 2
    ## [1] 3
    ## [1] 2

``` r
# make post drug a separate placebo condition
mergedDf$Drug2<-as.numeric(mergedDf$Drug)
mergedDf$Drug2[mergedDf$PostMDMASession==1]=3
# placebo pre-drug is 1
# placebo post-drug is 3
# drug is 2

mergedDf$Drug2<-as.factor(mergedDf$Drug2)
mergedDf <- within(mergedDf, Drug2 <- relevel(Drug2, ref = 1))
# set MDMA visits to post MDMA
mergedDf$PostMDMASession[mergedDf$Drug==1]=1

model <- lme(TDProp1 ~ MeanFD +Drug2+ RemTRs, random = ~ 1 | Subjects, data = mergedDf)
```

``` r
#meltedDf <- melt(mergedDf, id.vars = c("Subjects", "Task", "Dosage", "Session", "SpikesPercent", "MeanFD", "RemTRs", "Drug"))
#meltedDf <- within(meltedDf, variable <- relevel(variable, ref = 3))
#model <- lme(value ~ MeanFD + Drug + RemTRs+variable+variable*Drug+MeanFD*variable+RemTRs*variable, random = ~ 1 | Subjects/variable, data = meltedDf)
#summary(model)
#sjPlot::plot_model(model, type = "eff",show.values = T)
#
## great, now let's bootstrap it
## Set the number of bootstrap samples
#num_bootstrap_samples <- 
## initialize t vectors
#td1_d<-zeros(num_bootstrap_samples,1)
#td1_fd<-zeros(num_bootstrap_samples,1)
#td2_d<-zeros(num_bootstrap_samples,1)
#td2_fd<-zeros(num_bootstrap_samples,1)
#td3_d<-zeros(num_bootstrap_samples,1)
#td3_fd<-zeros(num_bootstrap_samples,1)
#td4_d<-zeros(num_bootstrap_samples,1)
#td4_fd<-zeros(num_bootstrap_samples,1)
#
## bootstrap loops
#set.seed(420)
#for (i in 1:num_bootstrap_samples){
#  # resample subjects instead of observations to be conserative
#   BootSubjs=sample(unique(mergedDf$Subjects),14,replace=T)
#   # Create an empty dataframe to store the resampled observations
#   bootSamp <- data.frame()
#   for (j in 1:length(BootSubjs)){
#       subject_obs <- meltedDf[meltedDf$Subjects == BootSubjs[j], ]
#       bootSamp <- rbind(bootSamp, subject_obs)
#   }
#  # fit model
#  model_i <- lme(value ~ MeanFD + Drug + RemTRs+variable+variable*Drug+variable*MeanFD, random = ~ 1 | Subjects/variable, data = bootSamp)
#  # get t-values
#  td1_d[i]=summary(model_i)$tTable[rownames(summary(model_i)$tTable) == "Drug1:variableTDProp1", "t-value"]
#  td1_fd[i]=summary(model_i)$tTable[rownames(summary(model_i)$tTable) == "MeanFD:variableTDProp1", "t-value"]
#  td3_d[i]=summary(model_i)$tTable[rownames(summary(model_i)$tTable) == "Drug1:variableTDProp3", "t-value"]
#  td3_fd[i]=summary(model_i)$tTable[rownames(summary(model_i)$tTable) == "MeanFD:variableTDProp3", "t-value"]
#  td4_d[i]=summary(model_i)$tTable[rownames(summary(model_i)$tTable) == "Drug1:variableTDProp4", "t-value"]
#  td4_fd[i]=summary(model_i)$tTable[rownames(summary(model_i)$tTable) == "MeanFD:variableTDProp4", "t-value"]
#}
#
#bootstrap_results=data.frame(cbind(td1_d,td1_fd,td3_d,td3_fd,td4_d,td4_fd))
#colnames(bootstrap_results)=c('td1_drug','td1_fd','td3_drug','td3_fd','td4_drug','td4_fd')
#plotdf=melt(bootstrap_results)
#
## Create box-and-whisker plots
#ggplot(plotdf, aes(x = variable, y = value,)) +
#  geom_boxplot(position = position_dodge(width = 0.8)) +
#  labs(x = "Model", y = "T-Values", title = "Bootstrap T-Values for Drug and MeanFD") +
#  theme_minimal(base_size = 16)+scale_fill_manual(values=c('blue','red'))
```

``` r
# autocorrelation
rs1=read.csv('~/Downloads/rs1_TAutoCorMerged.csv',header=F)
rs2=read.csv('~/Downloads/rs2_TAutoCorMerged.csv',header=F)
emo=read.csv('~/Downloads/emotion_TAutoCorMerged.csv',header=F)
gambling=read.csv('~/Downloads/gambling_TAutoCorMerged.csv',header=F)
wm=read.csv('~/Downloads/wm_TAutoCorMerged.csv',header=F)

rs1$Task='rs'
rs2$Task='rs2'
emo$Task='emotion'
gambling$Task='gambling'
wm$Task='wm'

# remove subj 4
rs1=rs1[-c(4),]
rs2=rs2[-c(4),]
emo=emo[-c(4),]
gambling=gambling[-c(4),]
wm=wm[-c(4),]

# this is going to be ugly but simple
# manually pair columns as sep. observations of baseline, placebo, 80, 120mg
rs1bv=data.frame(cbind(rs1$V1,rs1$V5))
rs1p=data.frame(cbind(rs1$V2,rs1$V6))
rs1m1=data.frame(cbind(rs1$V3,rs1$V7))
rs1m2=data.frame(cbind(rs1$V4,rs1$V8))
colnames(rs1bv)=c('AutoCor','RemTRs')
colnames(rs1p)=c('AutoCor','RemTRs')
colnames(rs1m1)=c('AutoCor','RemTRs')
colnames(rs1m2)=c('AutoCor','RemTRs')

rs2bv=data.frame(cbind(rs2$V1,rs2$V5))
rs2p=data.frame(cbind(rs2$V2,rs2$V6))
rs2m1=data.frame(cbind(rs2$V3,rs2$V7))
rs2m2=data.frame(cbind(rs2$V4,rs2$V8))
colnames(rs2bv)=c('AutoCor','RemTRs')
colnames(rs2p)=c('AutoCor','RemTRs')
colnames(rs2m1)=c('AutoCor','RemTRs')
colnames(rs2m2)=c('AutoCor','RemTRs')

emobv=data.frame(cbind(emo$V1,emo$V5))
emop=data.frame(cbind(emo$V2,emo$V6))
emom1=data.frame(cbind(emo$V3,emo$V7))
emom2=data.frame(cbind(emo$V4,emo$V8))
colnames(emobv)=c('AutoCor','RemTRs')
colnames(emop)=c('AutoCor','RemTRs')
colnames(emom1)=c('AutoCor','RemTRs')
colnames(emom2)=c('AutoCor','RemTRs')

gamblingbv=data.frame(cbind(gambling$V1,gambling$V5))
gamblingp=data.frame(cbind(gambling$V2,gambling$V6))
gamblingm1=data.frame(cbind(gambling$V3,gambling$V7))
gamblingm2=data.frame(cbind(gambling$V4,gambling$V8))
colnames(gamblingbv)=c('AutoCor','RemTRs')
colnames(gamblingp)=c('AutoCor','RemTRs')
colnames(gamblingm1)=c('AutoCor','RemTRs')
colnames(gamblingm2)=c('AutoCor','RemTRs')

wmbv=data.frame(cbind(wm$V1,wm$V5))
wmp=data.frame(cbind(wm$V2,wm$V6))
wmm1=data.frame(cbind(wm$V3,wm$V7))
wmm2=data.frame(cbind(wm$V4,wm$V8))
colnames(wmbv)=c('AutoCor','RemTRs')
colnames(wmp)=c('AutoCor','RemTRs')
colnames(wmm1)=c('AutoCor','RemTRs')
colnames(wmm2)=c('AutoCor','RemTRs')

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

emobv$Subjects=alff$SubjID
emop$Subjects=alff$SubjID
emom1$Subjects=alff$SubjID
emom2$Subjects=alff$SubjID

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

emobv$Task='emotion'
emop$Task='emotion'
emom1$Task='emotion'
emom2$Task='emotion'

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

emobv$Dosage='baseline'
emop$Dosage='Placebo'
emom1$Dosage='80mg'
emom2$Dosage='120mg'

gamblingbv$Dosage='baseline'
gamblingp$Dosage='Placebo'
gamblingm1$Dosage='80mg'
gamblingm2$Dosage='120mg'

wmbv$Dosage='baseline'
wmp$Dosage='Placebo'
wmm1$Dosage='80mg'
wmm2$Dosage='120mg'

# comibine all
allScans=rbind(rs1bv,rs1p,rs1m1,rs1m2,rs2bv,rs2p,rs2m1,rs2m2,emobv,emop,emom1,emom2,gamblingbv,gamblingp,gamblingm1,gamblingm2,wmbv,wmp,wmm1,wmm2)

# read in motion
mot=read.csv('~/Desktop/MDMA_spikes_summary.csv')

# motion merge
mergedDf=merge(mot,allScans,by=c('Subjects','Task','Dosage'))
#mergedDf=mergedDf[mergedDf$Dosage!='baseline',]
mergedDf$Drug=0
mergedDf$Drug[mergedDf$Dosage=="120mg"]=1
mergedDf$Drug[mergedDf$Dosage=="80mg"]=1
mergedDf$Drug=as.factor(mergedDf$Drug)
```

``` r
# combine complexity and props and autocor
mergedDfPropsComplAutoC=merge(mergedDfProps,mergedDf,by=c("Subjects","Task","Dosage","Session","MeanFD","SpikesPercent","RemTRs","Drug"))
```

``` r
# DMN seg
rs1=read.csv('~/Downloads/rs1_DMNSegMerged.csv',header=F)
rs2=read.csv('~/Downloads/rs2_DMNSegMerged.csv',header=F)
emo=read.csv('~/Downloads/emotion_DMNSegMerged.csv',header=F)
gambling=read.csv('~/Downloads/gambling_DMNSegMerged.csv',header=F)
wm=read.csv('~/Downloads/wm_DMNSegMerged.csv',header=F)
# set colnames
#colnames(rs1)=c('bvProp','pProp','m1Prop','m2Prop','bvTRs','pTRs','m1TRs','m2TRs')
#colnames(rs2)=c('bvProp','pProp','m1Prop','m2Prop','bvTRs','pTRs','m1TRs','m2TRs')
#colnames(emo)=c('bvProp','pProp','m1Prop','m2Prop','bvTRs','pTRs','m1TRs','m2TRs')
#colnames(gambling)=c('bvProp','pProp','m1Prop','m2Prop','bvTRs','pTRs','m1TRs','m2TRs')
#colnames(wm)=c('bvProp','pProp','m1Prop','m2Prop','bvTRs','pTRs','m1TRs','m2TRs')
rs1$Task='rs'
rs2$Task='rs2'
emo$Task='emotion'
gambling$Task='gambling'
wm$Task='wm'

# remove subj 4
rs1=rs1[-c(4),]
rs2=rs2[-c(4),]
emo=emo[-c(4),]
gambling=gambling[-c(4),]
wm=wm[-c(4),]

# this is going to be ugly but simple
# manually pair columns as sep. observations of baseline, placebo, 80, 120mg
rs1bv=data.frame(cbind(rs1$V1,rs1$V5))
rs1p=data.frame(cbind(rs1$V2,rs1$V6))
rs1m1=data.frame(cbind(rs1$V3,rs1$V7))
rs1m2=data.frame(cbind(rs1$V4,rs1$V8))
colnames(rs1bv)=c('DMNFC','RemTRs')
colnames(rs1p)=c('DMNFC','RemTRs')
colnames(rs1m1)=c('DMNFC','RemTRs')
colnames(rs1m2)=c('DMNFC','RemTRs')

rs2bv=data.frame(cbind(rs2$V1,rs2$V5))
rs2p=data.frame(cbind(rs2$V2,rs2$V6))
rs2m1=data.frame(cbind(rs2$V3,rs2$V7))
rs2m2=data.frame(cbind(rs2$V4,rs2$V8))
colnames(rs2bv)=c('DMNFC','RemTRs')
colnames(rs2p)=c('DMNFC','RemTRs')
colnames(rs2m1)=c('DMNFC','RemTRs')
colnames(rs2m2)=c('DMNFC','RemTRs')

emobv=data.frame(cbind(emo$V1,emo$V5))
emop=data.frame(cbind(emo$V2,emo$V6))
emom1=data.frame(cbind(emo$V3,emo$V7))
emom2=data.frame(cbind(emo$V4,emo$V8))
colnames(emobv)=c('DMNFC','RemTRs')
colnames(emop)=c('DMNFC','RemTRs')
colnames(emom1)=c('DMNFC','RemTRs')
colnames(emom2)=c('DMNFC','RemTRs')

gamblingbv=data.frame(cbind(gambling$V1,gambling$V5))
gamblingp=data.frame(cbind(gambling$V2,gambling$V6))
gamblingm1=data.frame(cbind(gambling$V3,gambling$V7))
gamblingm2=data.frame(cbind(gambling$V4,gambling$V8))
colnames(gamblingbv)=c('DMNFC','RemTRs')
colnames(gamblingp)=c('DMNFC','RemTRs')
colnames(gamblingm1)=c('DMNFC','RemTRs')
colnames(gamblingm2)=c('DMNFC','RemTRs')

wmbv=data.frame(cbind(wm$V1,wm$V5))
wmp=data.frame(cbind(wm$V2,wm$V6))
wmm1=data.frame(cbind(wm$V3,wm$V7))
wmm2=data.frame(cbind(wm$V4,wm$V8))
colnames(wmbv)=c('DMNFC','RemTRs')
colnames(wmp)=c('DMNFC','RemTRs')
colnames(wmm1)=c('DMNFC','RemTRs')
colnames(wmm2)=c('DMNFC','RemTRs')

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

emobv$Subjects=alff$SubjID
emop$Subjects=alff$SubjID
emom1$Subjects=alff$SubjID
emom2$Subjects=alff$SubjID

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

emobv$Task='emotion'
emop$Task='emotion'
emom1$Task='emotion'
emom2$Task='emotion'

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

emobv$Dosage='baseline'
emop$Dosage='Placebo'
emom1$Dosage='80mg'
emom2$Dosage='120mg'

gamblingbv$Dosage='baseline'
gamblingp$Dosage='Placebo'
gamblingm1$Dosage='80mg'
gamblingm2$Dosage='120mg'

wmbv$Dosage='baseline'
wmp$Dosage='Placebo'
wmm1$Dosage='80mg'
wmm2$Dosage='120mg'

# comibine all
allScans=rbind(rs1bv,rs1p,rs1m1,rs1m2,rs2bv,rs2p,rs2m1,rs2m2,emobv,emop,emom1,emom2,gamblingbv,gamblingp,gamblingm1,gamblingm2,wmbv,wmp,wmm1,wmm2)

# read in motion
mot=read.csv('~/Desktop/MDMA_spikes_summary.csv')

# motion merge
mergedDf=merge(mot,allScans,by=c('Subjects','Task','Dosage'))
#mergedDf=mergedDf[mergedDf$Dosage!='baseline',]
mergedDf$Drug=0
mergedDf$Drug[mergedDf$Dosage=="120mg"]=1
mergedDf$Drug[mergedDf$Dosage=="80mg"]=1
mergedDf$Drug=as.factor(mergedDf$Drug)
```

``` r
# combine complexity and props and autocor
mergedDfPropsComplAutoCdmnFC=merge(mergedDfPropsComplAutoC,mergedDf,by=c("Subjects","Task","Dosage","Session","MeanFD","SpikesPercent","RemTRs","Drug"))
```

``` r
# DMN mag
rs1=read.csv('~/Downloads/rs1_DMNMagMerged(4).csv',header=F)
rs2=read.csv('~/Downloads/rs2_DMNMagMerged(4).csv',header=F)
emo=read.csv('~/Downloads/emotion_DMNMagMerged(4).csv',header=F)
gambling=read.csv('~/Downloads/gambling_DMNMagMerged(4).csv',header=F)
wm=read.csv('~/Downloads/wm_DMNMagMerged(4).csv',header=F)
# set colnames
#colnames(rs1)=c('bvProp','pProp','m1Prop','m2Prop','bvTRs','pTRs','m1TRs','m2TRs')
#colnames(rs2)=c('bvProp','pProp','m1Prop','m2Prop','bvTRs','pTRs','m1TRs','m2TRs')
#colnames(emo)=c('bvProp','pProp','m1Prop','m2Prop','bvTRs','pTRs','m1TRs','m2TRs')
#colnames(gambling)=c('bvProp','pProp','m1Prop','m2Prop','bvTRs','pTRs','m1TRs','m2TRs')
#colnames(wm)=c('bvProp','pProp','m1Prop','m2Prop','bvTRs','pTRs','m1TRs','m2TRs')
rs1$Task='rs'
rs2$Task='rs2'
emo$Task='emotion'
gambling$Task='gambling'
wm$Task='wm'

# remove subj 4
rs1=rs1[-c(4),]
rs2=rs2[-c(4),]
emo=emo[-c(4),]
gambling=gambling[-c(4),]
wm=wm[-c(4),]

# this is going to be ugly but simple
# manually pair columns as sep. observations of baseline, placebo, 80, 120mg
rs1bv=data.frame(cbind(rs1$V1,rs1$V5))
rs1p=data.frame(cbind(rs1$V2,rs1$V6))
rs1m1=data.frame(cbind(rs1$V3,rs1$V7))
rs1m2=data.frame(cbind(rs1$V4,rs1$V8))
colnames(rs1bv)=c('DMNMag','RemTRs')
colnames(rs1p)=c('DMNMag','RemTRs')
colnames(rs1m1)=c('DMNMag','RemTRs')
colnames(rs1m2)=c('DMNMag','RemTRs')

rs2bv=data.frame(cbind(rs2$V1,rs2$V5))
rs2p=data.frame(cbind(rs2$V2,rs2$V6))
rs2m1=data.frame(cbind(rs2$V3,rs2$V7))
rs2m2=data.frame(cbind(rs2$V4,rs2$V8))
colnames(rs2bv)=c('DMNMag','RemTRs')
colnames(rs2p)=c('DMNMag','RemTRs')
colnames(rs2m1)=c('DMNMag','RemTRs')
colnames(rs2m2)=c('DMNMag','RemTRs')

emobv=data.frame(cbind(emo$V1,emo$V5))
emop=data.frame(cbind(emo$V2,emo$V6))
emom1=data.frame(cbind(emo$V3,emo$V7))
emom2=data.frame(cbind(emo$V4,emo$V8))
colnames(emobv)=c('DMNMag','RemTRs')
colnames(emop)=c('DMNMag','RemTRs')
colnames(emom1)=c('DMNMag','RemTRs')
colnames(emom2)=c('DMNMag','RemTRs')

gamblingbv=data.frame(cbind(gambling$V1,gambling$V5))
gamblingp=data.frame(cbind(gambling$V2,gambling$V6))
gamblingm1=data.frame(cbind(gambling$V3,gambling$V7))
gamblingm2=data.frame(cbind(gambling$V4,gambling$V8))
colnames(gamblingbv)=c('DMNMag','RemTRs')
colnames(gamblingp)=c('DMNMag','RemTRs')
colnames(gamblingm1)=c('DMNMag','RemTRs')
colnames(gamblingm2)=c('DMNMag','RemTRs')

wmbv=data.frame(cbind(wm$V1,wm$V5))
wmp=data.frame(cbind(wm$V2,wm$V6))
wmm1=data.frame(cbind(wm$V3,wm$V7))
wmm2=data.frame(cbind(wm$V4,wm$V8))
colnames(wmbv)=c('DMNMag','RemTRs')
colnames(wmp)=c('DMNMag','RemTRs')
colnames(wmm1)=c('DMNMag','RemTRs')
colnames(wmm2)=c('DMNMag','RemTRs')

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

emobv$Subjects=alff$SubjID
emop$Subjects=alff$SubjID
emom1$Subjects=alff$SubjID
emom2$Subjects=alff$SubjID

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

emobv$Task='emotion'
emop$Task='emotion'
emom1$Task='emotion'
emom2$Task='emotion'

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

emobv$Dosage='baseline'
emop$Dosage='Placebo'
emom1$Dosage='80mg'
emom2$Dosage='120mg'

gamblingbv$Dosage='baseline'
gamblingp$Dosage='Placebo'
gamblingm1$Dosage='80mg'
gamblingm2$Dosage='120mg'

wmbv$Dosage='baseline'
wmp$Dosage='Placebo'
wmm1$Dosage='80mg'
wmm2$Dosage='120mg'

# comibine all
allScans=rbind(rs1bv,rs1p,rs1m1,rs1m2,rs2bv,rs2p,rs2m1,rs2m2,emobv,emop,emom1,emom2,gamblingbv,gamblingp,gamblingm1,gamblingm2,wmbv,wmp,wmm1,wmm2)

# read in motion
mot=read.csv('~/Desktop/MDMA_spikes_summary.csv')

# motion merge
mergedDf=merge(mot,allScans,by=c('Subjects','Task','Dosage'))
mergedDfBVincl=mergedDf
mergedDf=mergedDf[mergedDf$Dosage!='baseline',]
mergedDf$Drug=0
mergedDf$Drug[mergedDf$Dosage=="120mg"]=1
mergedDf$Drug[mergedDf$Dosage=="80mg"]=1
mergedDf$Drug=as.factor(mergedDf$Drug)

mergedDfBVincl$Drug=0
mergedDfBVincl$Drug[mergedDfBVincl$Dosage=="120mg"]=1
mergedDfBVincl$Drug[mergedDfBVincl$Dosage=="80mg"]=1
mergedDfBVincl$Drug=as.factor(mergedDfBVincl$Drug)
```

``` r
# combine complexity and props and autocor
mergedDfPropsComplAutoCdmnFCdmnMag=merge(mergedDfPropsComplAutoCdmnFC,mergedDf,by=c("Subjects","Task","Dosage","Session","MeanFD","SpikesPercent","RemTRs","Drug"))
# version with baseline included
mergedDfPropsComplAutoCdmnFCdmnMagbv=merge(mergedDfPropsComplAutoCdmnFC,mergedDfBVincl,by=c("Subjects","Task","Dosage","Session","MeanFD","SpikesPercent","RemTRs","Drug"))

# remove data that needs to be removed (subs 6 and 10, <250 TRs)
mergedDfPropsComplAutoCdmnFCdmnMag=mergedDfPropsComplAutoCdmnFCdmnMag[mergedDfPropsComplAutoCdmnFCdmnMag$Subjects!='sub-MDMA006',]
mergedDfPropsComplAutoCdmnFCdmnMag=mergedDfPropsComplAutoCdmnFCdmnMag[mergedDfPropsComplAutoCdmnFCdmnMag$Subjects!='sub-MDMA010',]
mergedDfPropsComplAutoCdmnFCdmnMag=mergedDfPropsComplAutoCdmnFCdmnMag[mergedDfPropsComplAutoCdmnFCdmnMag$RemTRs>250,]

# change rs2 to rs for accurate task-modeling
mergedDfPropsComplAutoCdmnFCdmnMag$Task[mergedDfPropsComplAutoCdmnFCdmnMag$Task=='rs2']='rs'
mergedDfPropsComplAutoCdmnFCdmnMag$Task=as.factor(mergedDfPropsComplAutoCdmnFCdmnMag$Task)

# set rs to reference level
mergedDfPropsComplAutoCdmnFCdmnMag <- within(mergedDfPropsComplAutoCdmnFCdmnMag, Task <- relevel(Task, ref = 2))

### baseline included
# remove data that needs to be removed (subs 6 and 10, <250 TRs)
mergedDfPropsComplAutoCdmnFCdmnMagbv=mergedDfPropsComplAutoCdmnFCdmnMagbv[mergedDfPropsComplAutoCdmnFCdmnMagbv$Subjects!='sub-MDMA006',]
mergedDfPropsComplAutoCdmnFCdmnMagbv=mergedDfPropsComplAutoCdmnFCdmnMagbv[mergedDfPropsComplAutoCdmnFCdmnMagbv$Subjects!='sub-MDMA010',]
mergedDfPropsComplAutoCdmnFCdmnMagbv=mergedDfPropsComplAutoCdmnFCdmnMagbv[mergedDfPropsComplAutoCdmnFCdmnMagbv$RemTRs>250,]

# change rs2 to rs for accurate task-modeling
mergedDfPropsComplAutoCdmnFCdmnMagbv$Task[mergedDfPropsComplAutoCdmnFCdmnMagbv$Task=='rs2']='rs'
mergedDfPropsComplAutoCdmnFCdmnMagbv$Task=as.factor(mergedDfPropsComplAutoCdmnFCdmnMagbv$Task)

# set rs to reference level
mergedDfPropsComplAutoCdmnFCdmnMagbv <- within(mergedDfPropsComplAutoCdmnFCdmnMagbv, Task <- relevel(Task, ref = 2))
```

``` r
# model
td_model <- lme(TDProp1 ~ MeanFD + Drug+RemTRs+Task, random = ~ 1 | Subjects, data = mergedDfPropsComplAutoCdmnFCdmnMag)
ta_model <- lme(AutoCor ~ MeanFD + Drug+RemTRs+Task, random = ~ 1 | Subjects, data = mergedDfPropsComplAutoCdmnFCdmnMag)
dmnseg_model <- lme(DMNFC ~ MeanFD + Drug+RemTRs+Task, random = ~ 1 | Subjects, data = mergedDfPropsComplAutoCdmnFCdmnMag)
dmnmag_model <- lme(DMNMag ~ MeanFD + Drug+RemTRs+Task, random = ~ 1 | Subjects, data = mergedDfPropsComplAutoCdmnFCdmnMag)

# inclusive model
td_model <- lme(TDProp1 ~ MeanFD + Drug+RemTRs+Task+AutoCor+DMNFC, random = ~ 1 | Subjects, data = mergedDfPropsComplAutoCdmnFCdmnMag)
dmnmag_model <- lme(DMNMag ~ MeanFD + Drug+RemTRs+Task+AutoCor+DMNFC, random = ~ 1 | Subjects, data = mergedDfPropsComplAutoCdmnFCdmnMag)

# inclusive model, no drug
td_model <- lme(TDProp1 ~ MeanFD + Drug+RemTRs+Task+AutoCor+DMNFC, random = ~ 1 | Subjects, data = mergedDfPropsComplAutoCdmnFCdmnMagbv)
dmnmag_model <- lme(DMNMag ~ MeanFD + Drug+RemTRs+Task+AutoCor+DMNFC, random = ~ 1 | Subjects, data = mergedDfPropsComplAutoCdmnFCdmnMagbv)


# for correlation matrix
forCormat=data.frame(mergedDfPropsComplAutoCdmnFCdmnMag$TDProp1,mergedDfPropsComplAutoCdmnFCdmnMag$AutoCor,mergedDfPropsComplAutoCdmnFCdmnMag$DMNFC,mergedDfPropsComplAutoCdmnFCdmnMag$MeanFD,mergedDfPropsComplAutoCdmnFCdmnMag$RemTRs,mergedDfPropsComplAutoCdmnFCdmnMag$DMNMag)

colnames(forCormat)=c('TDProp1','AutoCor','DMNFC','MeanFD','RemTRs','DMNMag')

library(corrplot)
```

    ## corrplot 0.92 loaded

``` r
source("http://www.sthda.com/upload/rquery_cormat.r")

corrmatrix=rquery.cormat(forCormat,type="full")
```

![](Stats_n_Viz_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

``` r
corrplot(as.matrix(corrmatrix$r),method='number')
```

![](Stats_n_Viz_files/figure-gfm/unnamed-chunk-14-2.png)<!-- -->

``` r
# need to establish which DMN metrics are associated above and beyond drug effect
# DMN FC
dmnseg_model <- lme(DMNFC ~ MeanFD + Drug+RemTRs+Task+TDProp1+AutoCor+DMNMag, random = ~ 1 | Subjects, data = mergedDfPropsComplAutoCdmnFCdmnMag)
# DMN tdprops
td_model <- lme(TDProp1 ~ MeanFD + Drug+RemTRs+Task+DMNFC+AutoCor+DMNMag, random = ~ 1 | Subjects, data = mergedDfPropsComplAutoCdmnFCdmnMag)
# DMN Mag
dmnmag_model <- lme(DMNMag ~ MeanFD + Drug+RemTRs+Task+DMNFC+AutoCor+TDProp1, random = ~ 1 | Subjects, data = mergedDfPropsComplAutoCdmnFCdmnMag)
# DMN autocor
dmnAC_model <- lme(AutoCor ~ MeanFD + Drug+RemTRs+Task+DMNFC+DMNMag+TDProp1, random = ~ 1 | Subjects, data = mergedDfPropsComplAutoCdmnFCdmnMag)

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
mergedDfPropsComplAutoCdmnFCdmnMag_notask=mergedDfPropsComplAutoCdmnFCdmnMag[mergedDfPropsComplAutoCdmnFCdmnMag$Task=='rs',]

# Fit logistic regression models
model1 <- glm(Drug ~ MeanFD + RemTRs + DMNFC+AutoCor, data = mergedDfPropsComplAutoCdmnFCdmnMag_notask, family = binomial)
model2 <- glm(Drug ~ MeanFD + RemTRs + DMNFC + AutoCor+DMNMag+TDProp1, data = mergedDfPropsComplAutoCdmnFCdmnMag_notask, family = binomial)

# Predict probabilities
prob1 <- predict(model1, type = "response")
prob2 <- predict(model2, type = "response")

# Create a combined data frame for all models
df <- data.frame(
  labels = as.numeric(rep(mergedDfPropsComplAutoCdmnFCdmnMag_notask$Drug, 2)),
  predictions = c(prob1, prob2),
  model = factor(rep(c("DMN Correlations", "+DMN Propagations"), each = nrow(mergedDfPropsComplAutoCdmnFCdmnMag_notask)))
)

# Generate the ROC plot
ggplot(df, aes(m = predictions, d = labels, color = model)) + 
  geom_roc(n.cuts = 0, labels = FALSE) + 
  ylim(0, 1) + ylab('True Positive Rate') +xlab('False Positive Rate')+
  ggtitle("ROC Curves for Classifying MDMA") + 
  theme_minimal(base_size=18) + 
  scale_color_manual(values = c("#09416b","#c12139"))+
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray")+
  theme(legend.position = "none")
```

    ## Warning in verify_d(data$d): D not labeled 0/1, assuming 1 = 0 and 2 = 1!

![](Stats_n_Viz_files/figure-gfm/unnamed-chunk-14-3.png)<!-- -->

``` r
# Calculate AUC for each model
roc1 <- roc(mergedDfPropsComplAutoCdmnFCdmnMag_notask$Drug, prob1)
```

    ## Setting levels: control = 0, case = 1

    ## Setting direction: controls < cases

``` r
roc2 <- roc(mergedDfPropsComplAutoCdmnFCdmnMag_notask$Drug, prob2)
```

    ## Setting levels: control = 0, case = 1
    ## Setting direction: controls < cases

``` r
# Print AUC values
auc1 <- auc(roc1)
auc2 <- auc(roc2)


print(paste("AUC for DMN Correlations:", auc1))
```

    ## [1] "AUC for DMN Correlations: 0.806481481481482"

``` r
print(paste("AUC for Full Model:", auc2))
```

    ## [1] "AUC for Full Model: 0.881481481481482"

``` r
# Calculate AUC difference between full and reduced models
auc_diff <- auc2 - auc1
```

``` r
# Temporary comment out: output is making rstudio wonky




# make equivalent AUC calculations on permuted data
# initialize AUC difference vectors
#auc_diffs <- rep(NA, 1000)
#
## 1. permute each DMN variable (and FD)
#set.seed(1)
#for (i in 1:1000){
#  print(i)
#  # permute DMNMag
#  mergedDfPropsComplAutoCdmnFCdmnMag_notask$DMNMag_perm <- sample(mergedDfPropsComplAutoCdmnFCdmnMag_notask$DMNMag)
#  # permute TDProp1
#  mergedDfPropsComplAutoCdmnFCdmnMag_notask$TDProp1_perm <- sample(mergedDfPropsComplAutoCdmnFCdmnMag_notask$TDProp1)
# 
#  # Fit logistic regression models
#  model1 <- glm(Drug ~ MeanFD + RemTRs + DMNFC+AutoCor, data = mergedDfPropsComplAutoCdmnFCdmnMag_notask, family = #binomial)
#  model2_perm <- glm(Drug ~ MeanFD + RemTRs + DMNFC+AutoCor+TDProp1_perm+DMNMag_perm, data = #mergedDfPropsComplAutoCdmnFCdmnMag_notask, family = binomial)
#   
#  # 3. calculate AUC difference between full and reduced models with permuted data
#  roc1 <- roc(mergedDfPropsComplAutoCdmnFCdmnMag_notask$Drug, predict(model1, type = "response"))
#  roc2_perm <- roc(mergedDfPropsComplAutoCdmnFCdmnMag_notask$Drug, predict(model2_perm, type = "response"))
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
# great, now let's bootstrap them
# Set the number of bootstrap samples
num_bootstrap_samples <- 1000
# initialize t vectors
td_d<-rep(0,num_bootstrap_samples)
td_fd<-rep(0,num_bootstrap_samples)
ta_d<-rep(0,num_bootstrap_samples)
ta_fd<-rep(0,num_bootstrap_samples)
ds_d<-rep(0,num_bootstrap_samples)
ds_fd<-rep(0,num_bootstrap_samples)
dm_d<-rep(0,num_bootstrap_samples)
dm_fd<-rep(0,num_bootstrap_samples)

# bootstrap loops
set.seed(1)
for (i in 1:num_bootstrap_samples){
  # resample data
  data=mergedDfPropsComplAutoCdmnFCdmnMag[sample(nrow(mergedDfPropsComplAutoCdmnFCdmnMag), replace = TRUE), ]
  # fit on all models
  td_model <- lme(TDProp1 ~ MeanFD + Drug+RemTRs+Task, random = ~ 1 | Subjects, data = data)
  ta_model <- lme(AutoCor ~ MeanFD + Drug+RemTRs+Task, random = ~ 1 | Subjects, data = data)
  ds_model <- lme(DMNFC ~ MeanFD + Drug+RemTRs+Task, random = ~ 1 | Subjects, data = data)
  dm_model <- lme(DMNMag ~ MeanFD + Drug+RemTRs+Task, random = ~ 1 | Subjects, data = data)
  # get t-values
  td_d[i]=summary(td_model)$tTable[rownames(summary(td_model)$tTable) == "Drug1", "t-value"]
  td_fd[i]=summary(td_model)$tTable[rownames(summary(td_model)$tTable) == "MeanFD", "t-value"]
  ta_d[i]=summary(ta_model)$tTable[rownames(summary(ta_model)$tTable) == "Drug1", "t-value"]
  ta_fd[i]=summary(ta_model)$tTable[rownames(summary(ta_model)$tTable) == "MeanFD", "t-value"]
  ds_d[i]=summary(ds_model)$tTable[rownames(summary(ds_model)$tTable) == "Drug1", "t-value"]
  ds_fd[i]=summary(ds_model)$tTable[rownames(summary(ds_model)$tTable) == "MeanFD", "t-value"]
  dm_d[i]=summary(dm_model)$tTable[rownames(summary(dm_model)$tTable) == "Drug1", "t-value"]
  dm_fd[i]=summary(dm_model)$tTable[rownames(summary(dm_model)$tTable) == "MeanFD", "t-value"]
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
  summarize(avg_t_value = mean(abs(tstat), na.rm = TRUE))

# Reorder the Model factor based on the average t-values
bootstrap_results_FD$Model <- factor(bootstrap_results_FD$Model, 
                                     levels = average_t_values$Model[order(average_t_values$avg_t_value)])


# Calculate the average t-value for each Model category
average_t_values <- bootstrap_results_Drug %>%
  group_by(Model) %>%
  summarize(avg_t_value = mean(abs(tstat), na.rm = TRUE))

# Reorder the Model factor based on the average t-values
bootstrap_results_Drug$Model <- factor(bootstrap_results_Drug$Model, 
                                     levels = average_t_values$Model[order(average_t_values$avg_t_value)])

library(ggdist)
# Generate the plot
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

![](Stats_n_Viz_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

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
    labs(x = "Model", y = "T-Values", title = "Bootstrap T-Values for MDMA effect") +
    theme_minimal(base_size = 18) +
    geom_hline(yintercept = 0, linetype = "dashed")+
    # just to prevent extra x-axis expansion
    coord_cartesian(xlim = c(1, length(unique(bootstrap_results_Drug$Model))))+
    theme(legend.position = "none")+ylim(c(-10.5,10.5))
```

![](Stats_n_Viz_files/figure-gfm/unnamed-chunk-16-2.png)<!-- -->

``` r
# now compare with VAS
vas=read.csv('~/Downloads/MDMA_VAS_ASC_forAdam_03062023.csv')
# correct subject naming
vas$Subjects <- paste0("sub-MDMA", sprintf("%03d", vas$Subjects))
VASmerge=merge(vas,mergedDfPropsComplAutoCdmnFCdmnMag,by=c("Subjects","Dosage"))
```

``` r
###### spin tests
# load in y7 alignment left
y7align_L=read.csv('~/Downloads/perm_vs_obs_DMN_Y7_L.csv')
y7align_R=read.csv('~/Downloads/perm_vs_obs_DMN_Y7_R.csv')
# add a real vs. permuted value
y7align_L$Observed=0
y7align_R$Observed=0
y7align_L$Observed[10001]=1
y7align_R$Observed[10001]=1
# add right vs. left
y7align_L$Side='Left'
y7align_R$Side='Right'
y7plotdf=rbind(y7align_L,y7align_R)

# subset spun and nonspun df
spunL=y7align_L[1:10000,]
ObsL=y7align_L[10001,]
spunR=y7align_R[1:10000,]
ObsR=y7align_R[10001,]
  
# left hemi
ggplot(spunL,aes(x=Var1))+geom_density(size=1.5)+geom_vline(xintercept = ObsL$Var1,size=2,color='#BC3754')+theme_classic(base_size=23)+ylab('')+xlab('T-Statistics')+guides(y="none")+theme(axis.text = element_text(size=22))+ggtitle('DMN localization to Yeo7 Boundary: Left hemisphere')
```

    ## Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
    ## ℹ Please use `linewidth` instead.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

![](Stats_n_Viz_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

``` r
# right hemi
ggplot(spunR,aes(x=Var1))+geom_density(size=1.5)+geom_vline(xintercept = ObsR$Var1,size=2,color='#BC3754')+theme_classic(base_size=23)+ylab('')+xlab('T-Statistics')+guides(y="none")+theme(axis.text = element_text(size=22))+ggtitle('DMN localization to Yeo7 Boundary: Right hemisphere')
```

![](Stats_n_Viz_files/figure-gfm/unnamed-chunk-18-2.png)<!-- -->

``` r
# now for biol psych rois
bpalign_L=read.csv('~/Downloads/perm_vs_obs_DMN_BiolPsych_L.csv')
bpalign_R=read.csv('~/Downloads/perm_vs_obs_DMN_BiolPsych_R.csv')
# add a real vs. permuted value
bpalign_L$Observed=0
bpalign_R$Observed=0
bpalign_L$Observed[10001]=1
bpalign_R$Observed[10001]=1
# add right vs. left
bpalign_L$Side='Left'
bpalign_R$Side='Right'

# subset spun and nonspun df
spunL=bpalign_L[1:10000,]
ObsL=bpalign_L[10001,]
spunR=bpalign_R[1:10000,]
ObsR=bpalign_R[10001,]
  
# left hemi
ggplot(spunL,aes(x=Var1))+geom_density(size=1.5)+geom_vline(xintercept = ObsL$Var1,size=2,color='#BC3754')+theme_classic(base_size=23)+ylab('')+xlab('T-Statistics')+guides(y="none")+theme(axis.text = element_text(size=22))+ggtitle('DMN localization to Biol. Psych. ROIs: Left hemisphere')
```

![](Stats_n_Viz_files/figure-gfm/unnamed-chunk-18-3.png)<!-- -->

``` r
# right hemi
ggplot(spunR,aes(x=Var1))+geom_density(size=1.5)+geom_vline(xintercept = ObsR$Var1,size=2,color='#BC3754')+theme_classic(base_size=23)+ylab('')+xlab('T-Statistics')+guides(y="none")+theme(axis.text = element_text(size=22))+ggtitle('DMN localization to Biol. Psych. ROIs: Right hemisphere')
```

![](Stats_n_Viz_files/figure-gfm/unnamed-chunk-18-4.png)<!-- -->

``` r
# calculate change in DAS to investigate inter-psychedelic-session variability
# initialize change columns
DAS_Change_Subject=c('')
DAS_Change_Dosage=c('')
DAS_Change_ChangeTA=0
DAS_Change_ChangeFC=0
DAS_Change_ChangeBUP=0
DAS_Change_ChangeMag=0

# isolate drug sessions
VASmerge_DrugSeshs=VASmerge[VASmerge$Drug==1,]
VASmerge_Plac=VASmerge[VASmerge$Drug==0,]
eachSubj=unique(VASmerge$Subjects)

# create change from placebo column
for (i in 1:length(eachSubj)){
  # extract drug sessions
  SubjDrug=VASmerge_DrugSeshs[VASmerge_DrugSeshs$Subjects==eachSubj[i],]
  Subj80=SubjDrug[SubjDrug$Dosage=='80mg',]
  Subj120=SubjDrug[SubjDrug$Dosage=='120mg',]
  # extract placebo session
  SubjPlac=VASmerge_Plac[VASmerge_Plac$Subjects==eachSubj[i],]

  # calculate mean change across scans (TA)
  TA_Plac=mean(SubjPlac$AutoCor)
  TA_80=mean(Subj80$AutoCor)
  TA_120=mean(Subj120$AutoCor)
  SubjChange_TA_80=TA_Plac-TA_80
  SubjChange_TA_120=TA_Plac-TA_120
  
  # calculate mean change across scans (DMN integration)
  FC_Plac=mean(SubjPlac$DMNFC)
  FC_80=mean(Subj80$DMNFC)
  FC_120=mean(Subj120$DMNFC)
  SubjChange_FC_80=FC_Plac-FC_80
  SubjChange_FC_120=FC_Plac-FC_120
  
  # calculate mean change across scans (BUP)
  BUP_Plac=mean(SubjPlac$TDProp1)
  BUP_80=mean(Subj80$TDProp1)
  BUP_120=mean(Subj120$TDProp1)
  SubjChange_BUP_80=BUP_Plac-BUP_80
  SubjChange_BUP_120=BUP_Plac-BUP_120
  
  # calculate mean change across scans (Magnitude)
  Mag_Plac=mean(SubjPlac$DMNMag)
  Mag_80=mean(Subj80$DMNMag)
  Mag_120=mean(Subj120$DMNMag)
  SubjChange_Mag_80=Mag_Plac-Mag_80
  SubjChange_Mag_120=Mag_Plac-Mag_120
  
  # put into VAS_change, 80mg
  DAS_Change_Subject<-c(DAS_Change_Subject,eachSubj[i])
  DAS_Change_Dosage<-c(DAS_Change_Dosage,'80mg')
  DAS_Change_ChangeFC<-c(DAS_Change_ChangeFC,SubjChange_TA_80)
  DAS_Change_ChangeTA<-c(DAS_Change_ChangeTA,SubjChange_FC_80)
  DAS_Change_ChangeBUP<-c(DAS_Change_ChangeBUP,SubjChange_BUP_80)
  DAS_Change_ChangeMag<-c(DAS_Change_ChangeMag,SubjChange_Mag_80)
  # 120
  DAS_Change_Subject<-c(DAS_Change_Subject,eachSubj[i])
  DAS_Change_Dosage<-c(DAS_Change_Dosage,'120mg')
  DAS_Change_ChangeFC<-c(DAS_Change_ChangeFC,SubjChange_TA_120)
  DAS_Change_ChangeTA<-c(DAS_Change_ChangeTA,SubjChange_FC_120)
  DAS_Change_ChangeBUP<-c(DAS_Change_ChangeBUP,SubjChange_BUP_120)
  DAS_Change_ChangeMag<-c(DAS_Change_ChangeMag,SubjChange_Mag_120)
}

DAS_changeDF=data.frame(DAS_Change_Subject,DAS_Change_Dosage,DAS_Change_ChangeFC,DAS_Change_ChangeTA,DAS_Change_ChangeBUP,DAS_Change_ChangeMag)
# drop first row (initialization row)
DAS_changeDF=DAS_changeDF[-c(1),]
colnames(DAS_changeDF)<-c('Subjects','Dosage','FC_Decrease','TA_Decrease','BUP_Decrease','Mag_Decrease')
# merge with DAS scores
DAS_changeDF=merge(DAS_changeDF,vas,by=c("Subjects","Dosage"))
# omit NA rows
DAS_changeDF=DAS_changeDF[DAS_changeDF$FC_Decrease!='NaN',]



# initialize output correlation stats
corvec=NULL
pvec=NULL
colNameVec=NULL
counter=1
# plot change from placebo with reported DAS scores
for (i in 41:57){
  plot(DAS_changeDF$BUP_Decrease,DAS_changeDF[,i],main=colnames(DAS_changeDF)[i])
  # note das scores are not normally distriuted, even from drug sessions
  a=cor.test(DAS_changeDF$BUP_Decrease,DAS_changeDF[,i])
  corvec[counter]=a$estimate
  pvec[counter]=a$p.value
  colNameVec[counter]=colnames(DAS_changeDF)[i]
  counter=counter+1
}
```

![](Stats_n_Viz_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->![](Stats_n_Viz_files/figure-gfm/unnamed-chunk-19-2.png)<!-- -->![](Stats_n_Viz_files/figure-gfm/unnamed-chunk-19-3.png)<!-- -->![](Stats_n_Viz_files/figure-gfm/unnamed-chunk-19-4.png)<!-- -->![](Stats_n_Viz_files/figure-gfm/unnamed-chunk-19-5.png)<!-- -->![](Stats_n_Viz_files/figure-gfm/unnamed-chunk-19-6.png)<!-- -->![](Stats_n_Viz_files/figure-gfm/unnamed-chunk-19-7.png)<!-- -->![](Stats_n_Viz_files/figure-gfm/unnamed-chunk-19-8.png)<!-- -->![](Stats_n_Viz_files/figure-gfm/unnamed-chunk-19-9.png)<!-- -->![](Stats_n_Viz_files/figure-gfm/unnamed-chunk-19-10.png)<!-- -->![](Stats_n_Viz_files/figure-gfm/unnamed-chunk-19-11.png)<!-- -->![](Stats_n_Viz_files/figure-gfm/unnamed-chunk-19-12.png)<!-- -->![](Stats_n_Viz_files/figure-gfm/unnamed-chunk-19-13.png)<!-- -->![](Stats_n_Viz_files/figure-gfm/unnamed-chunk-19-14.png)<!-- -->![](Stats_n_Viz_files/figure-gfm/unnamed-chunk-19-15.png)<!-- -->![](Stats_n_Viz_files/figure-gfm/unnamed-chunk-19-16.png)<!-- -->![](Stats_n_Viz_files/figure-gfm/unnamed-chunk-19-17.png)<!-- -->

``` r
# set some colorations to be consistent in new df by subject
DAS_changeDF$Subjects<-as.factor(DAS_changeDF$Subjects)
DAS_changeDF$Colors <- generated_colors[DAS_changeDF$Subjects]
DAS_changeDF_people <- merge(DAS_changeDF, unique_pairs, by = 'Subjects', all.x = TRUE)
DAS_changeDF_people$SubjectsCols <- factor(DAS_changeDF_people$Subjects, levels = names(generated_colors))
```

``` r
### figure 3B - coarse
# make a vector of DAS scale names \n is a newline for the ggplots
dascNames <- c(
  'Oceaninc\nBoundlessness',
  'Dread of Ego\nDissolution',
  'Visionary\nRestructuralization',
  'Auditory\nAlterations',
  'Vigilance\nReduction'
)

dascNamesAcr <- c(
  'O.B.',
  'D.E.D.',
  'V.Res.',
  'A.A.',
  'V.Red.'
)


##### Propagation directions
# initialize correlation vector for dasc
BUP_corvec=NULL
BUP_pvec=NULL
colNameVec=NULL
# initialize counter
counter=1
for (i in 44:48){
  a=cor.test(DAS_changeDF$BUP_Decrease,DAS_changeDF[,i])
  BUP_corvec[counter] <- a$estimate
  BUP_pvec[counter] <- a$p.value
  counter=counter+1
}


#########
# equivalent for AutoCor
#########

TA_corvec=NULL
TA_pvec=NULL
# initialize counter
counter=1
for (i in 44:48){
  a=cor.test(DAS_changeDF$TA_Decrease,DAS_changeDF[,i])
  TA_corvec[counter] <- a$estimate
  TA_pvec[counter] <- a$p.value
  counter=counter+1
}


#########
# equivalent for DMN seg
#########

S_corvec=NULL
S_pvec=NULL
# initialize counter
counter=1
for (i in 44:48){
  a=cor.test(DAS_changeDF$FC_Decrease,DAS_changeDF[,i])
  S_corvec[counter] <- a$estimate
  S_pvec[counter] <- a$p.value
  counter=counter+1
}

#########
# equivalent for DMN Magnitude
#########

M_corvec=NULL
M_pvec=NULL
# initialize counter
counter=1
for (i in 44:48){
  a=cor.test(DAS_changeDF$Mag_Decrease,DAS_changeDF[,i])
  M_corvec[counter] <- a$estimate
  M_pvec[counter] <- a$p.value
  counter=counter+1
}

# correct all for multiple comparisons
allPs=c(S_pvec,TA_pvec,BUP_pvec,M_pvec)
allPs_MC=p.adjust(allPs,method='fdr')

# pull out S pvecs
S_pvec_mc=allPs_MC[1:5]
TA_pvec_mc=allPs_MC[6:10]
BUP_pvec_mc=allPs_MC[11:15]
M_pvec_mc=allPs_MC[16:20]

# BAR PLOTS FOR ALL
# DMN SEGREGATION
# Create a dataframe with the values and column names
df <- data.frame(Values = S_corvec, pvals=S_pvec)

# make a color vector by significance
colorvec=rep('NS',5)
colorvec[S_pvec<.05]='Uncr'
colorvec[S_pvec_mc<.05]='FDR'
df$Sig <- factor(colorvec, levels = c('NS', 'Uncr', 'FDR'))
# Specify color scale manually
colors <- c('NS' = '#002642', 'Uncr' = '#EF9500', 'FDR' = '#840032')
df$ColumnNamesAcr<-c(dascNames)
# linear version
ggplot(df, aes(x = Values, y = dascNames,fill=Sig)) +
  geom_bar(size=3,stat = "identity")+
  scale_fill_manual(values = colors) +
  labs(y = "", x = "Subscale Correlation",title="Change in Integration")+theme_minimal(base_size=17)+xlim(c(-.3,.8))+
  theme(legend.position='none')
```

![](Stats_n_Viz_files/figure-gfm/unnamed-chunk-20-1.png)<!-- -->

``` r
# DMN AUTOCOTR
# Create a dataframe with the values and column names
df <- data.frame(Values = TA_corvec, pvals=TA_pvec)

# make a color vector by significance
colorvec=rep('NS',5)
colorvec[TA_pvec<.05]='Uncr'
colorvec[TA_pvec_mc<.05]='FDR'
df$Sig <- factor(colorvec, levels = c('NS', 'Uncr', 'FDR'))
# Specify color scale manually
colors <- c('NS' = '#002642', 'Uncr' = '#EF9500', 'FDR' = '#840032')
df$ColumnNamesAcr<-c(dascNamesAcr)
# bar version
ggplot(df, aes(x = Values, y = ColumnNamesAcr,fill=Sig)) +
  geom_bar(size=3,stat = "identity")+
  scale_fill_manual(values = colors) +
  labs(y = "", x = "Subscale Correlation",title="Change in Autocorrelation")+theme_minimal(base_size=17)+xlim(c(-.3,.8))+
  theme(legend.position='none')
```

![](Stats_n_Viz_files/figure-gfm/unnamed-chunk-20-2.png)<!-- -->

``` r
# DMN BUP
# Create a dataframe with the values and column names
df <- data.frame(Values = BUP_corvec, pvals=BUP_pvec)

# make a color vector by significance
colorvec=rep('NS',5)
colorvec[BUP_pvec<.05]='Uncr'
colorvec[BUP_pvec_mc<.05]='FDR'
df$Sig <- factor(colorvec, levels = c('NS', 'Uncr', 'FDR'))
# Specify color scale manually
colors <- c('NS' = '#002642', 'Uncr' = '#EF9500', 'FDR' = '#840032')
df$ColumnNamesAcr<-c(dascNamesAcr)
# bar version
ggplot(df, aes(x = Values, y = ColumnNamesAcr,fill=Sig)) +
  geom_bar(size=3,stat = "identity")+
  scale_fill_manual(values = colors) +
  labs(y = "", x = "Subscale Correlation",title="Change in Prop. Direction")+theme_minimal(base_size=17)+xlim(c(-.3,.8))+
  theme(legend.position='none')
```

![](Stats_n_Viz_files/figure-gfm/unnamed-chunk-20-3.png)<!-- -->

``` r
# DMN Mag
# Create a dataframe with the values and column names
df <- data.frame(Values = M_corvec, pvals=M_pvec)

# make a color vector by significance
colorvec=rep('NS',5)
colorvec[M_pvec<.05]='Uncr'
colorvec[M_pvec_mc<.05]='FDR'
df$Sig <- factor(colorvec, levels = c('NS', 'Uncr', 'FDR'))
# Specify color scale manually
colors <- c('NS' = '#002642', 'Uncr' = '#EF9500', 'FDR' = '#840032')
df$ColumnNamesAcr<-c(dascNamesAcr)
# bar version
ggplot(df, aes(x = Values, y = ColumnNamesAcr,fill=Sig)) +
  geom_bar(size=3,stat = "identity")+
  scale_fill_manual(values = colors) +
  labs(y = "", x = "Subscale Correlation",title="Change in Prop. Magnitudes")+theme_minimal(base_size=17)+xlim(c(-.4,.9))+
  theme(legend.position='none')
```

![](Stats_n_Viz_files/figure-gfm/unnamed-chunk-20-4.png)<!-- -->

``` r
# pull out BUP~DED correlation
ggplot(data=DAS_changeDF_people,aes(y=dascscore_ded,x=BUP_Decrease))+
  geom_point(size=4,aes(color = People))+
  labs(x = "Decrease in Bottom-up %", y = "Increase in Dread of Ego Dissolution", color = "People") +
  scale_color_manual(values = generated_colors) +
  theme_minimal(base_size = 16)
```

![](Stats_n_Viz_files/figure-gfm/unnamed-chunk-20-5.png)<!-- -->

``` r
cor.test(DAS_changeDF_people$dascscore_ded,DAS_changeDF_people$BUP_Decrease)
```

    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  DAS_changeDF_people$dascscore_ded and DAS_changeDF_people$BUP_Decrease
    ## t = 4.5982, df = 20, p-value = 0.0001741
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  0.4231442 0.8742441
    ## sample estimates:
    ##       cor 
    ## 0.7168634

``` r
### figure 3B - granular
# make a vector of DAS scale names \n is a newline for the ggplots
dascNames=c('Experience of Unity','Spirtual Experience','Blissful State','Insightfulness','Disembodiment','Impaired Control','Anxiety','Complex Imagery','Elementary Imagery','Synsthesia','Changed meaning\n of precepts')

dascNamesAcr <- c(
  'E.U.',
  'S.E.',
  'B.S.',
  'Ins.',
  'Disemb.',
  'I.C.',
  'Anx.',
  'C.Im.',
  'E.Im.',
  'Syn.',
  'C.M.O.P.'
)
##### Propagation directions
# initialize correlation vector for dasc
BUP_corvec=NULL
BUP_pvec=NULL
colNameVec=NULL
# initialize counter
counter=1
for (i in 50:60){
  a=cor.test(DAS_changeDF$BUP_Decrease,DAS_changeDF[,i])
  BUP_corvec[counter] <- a$estimate
  BUP_pvec[counter] <- a$p.value
  counter=counter+1
}


#########
# equivalent for AutoCor
#########

TA_corvec=NULL
TA_pvec=NULL
# initialize counter
counter=1
for (i in 50:60){
  a=cor.test(DAS_changeDF$TA_Decrease,DAS_changeDF[,i])
  TA_corvec[counter] <- a$estimate
  TA_pvec[counter] <- a$p.value
  counter=counter+1
}


#########
# equivalent for DMN seg
#########

S_corvec=NULL
S_pvec=NULL
# initialize counter
counter=1
for (i in 50:60){
  a=cor.test(DAS_changeDF$FC_Decrease,DAS_changeDF[,i])
  S_corvec[counter] <- a$estimate
  S_pvec[counter] <- a$p.value
  counter=counter+1
}

#########
# equivalent for DMN Magnitude
#########

M_corvec=NULL
M_pvec=NULL
# initialize counter
counter=1
for (i in 50:60){
  a=cor.test(DAS_changeDF$Mag_Decrease,DAS_changeDF[,i])
  M_corvec[counter] <- a$estimate
  M_pvec[counter] <- a$p.value
  counter=counter+1
}

# correct all for multiple comparisons
allPs=c(S_pvec,TA_pvec,BUP_pvec,M_pvec)
allPs_MC=p.adjust(allPs,method='fdr')

# pull out S pvecs
S_pvec_mc=allPs_MC[1:11]
TA_pvec_mc=allPs_MC[12:22]
BUP_pvec_mc=allPs_MC[23:33]
M_pvec_mc=allPs_MC[34:44]

# BAR PLOTS FOR ALL
# DMN SEGREGATION
# Create a dataframe with the values and column names
df <- data.frame(Values = S_corvec, pvals=S_pvec)

# make a color vector by significance
colorvec=rep('NS',11)
colorvec[S_pvec<.05]='Uncr'
colorvec[S_pvec_mc<.05]='FDR'
df$Sig <- factor(colorvec, levels = c('NS', 'Uncr', 'FDR'))
# Specify color scale manually
colors <- c('NS' = '#002642', 'Uncr' = '#EF9500', 'FDR' = '#840032')
df$ColumnNamesAcr<-c(dascNames)
# linear version
ggplot(df, aes(x = Values, y = dascNames,fill=Sig)) +
  geom_bar(size=3,stat = "identity")+
  scale_fill_manual(values = colors) +
  labs(y = "", x = "Subscale Correlation",title="Change in Integration")+theme_minimal(base_size=17)+xlim(c(-.3,.8))+
  theme(legend.position='none')
```

![](Stats_n_Viz_files/figure-gfm/unnamed-chunk-21-1.png)<!-- -->

``` r
# DMN AUTOCOR
# Create a dataframe with the values and column names
df <- data.frame(Values = TA_corvec, pvals=TA_pvec)

# make a color vector by significance
colorvec=rep('NS',11)
colorvec[TA_pvec<.05]='Uncr'
colorvec[TA_pvec_mc<.05]='FDR'
df$Sig <- factor(colorvec, levels = c('NS', 'Uncr', 'FDR'))
# Specify color scale manually
colors <- c('NS' = '#002642', 'Uncr' = '#EF9500', 'FDR' = '#840032')
# extra step needs to be done to maintain same row ordering as previous plot
df$ColumnNamesAcr<-c(dascNamesAcr)
df$ColumnNamesAcr<-factor(df$ColumnNamesAcr,levels=c('Anx.','B.S.','C.M.O.P.','C.Im.','Disemb.','E.Im.','E.U.','I.C.','Ins.','S.E.','Syn.'))
# bar version
ggplot(df, aes(x = Values, y = ColumnNamesAcr,fill=Sig)) +
  geom_bar(size=3,stat = "identity")+
  scale_fill_manual(values = colors) +
  labs(y = "", x = "Subscale Correlation",title="Change in Autocorrelation")+theme_minimal(base_size=17)+xlim(c(-.3,.8))+
  theme(legend.position='none')
```

![](Stats_n_Viz_files/figure-gfm/unnamed-chunk-21-2.png)<!-- -->

``` r
# DMN BUP
# Create a dataframe with the values and column names
df <- data.frame(Values = BUP_corvec, pvals=BUP_pvec)

# make a color vector by significance
colorvec=rep('NS',11)
colorvec[BUP_pvec<.05]='Uncr'
colorvec[BUP_pvec_mc<.05]='FDR'
df$Sig <- factor(colorvec, levels = c('NS', 'Uncr', 'FDR'))
# Specify color scale manually
colors <- c('NS' = '#002642', 'Uncr' = '#EF9500', 'FDR' = '#840032')
# extra step needs to be done to maintain same row ordering as previous plot
df$ColumnNamesAcr<-c(dascNamesAcr)
df$ColumnNamesAcr<-factor(df$ColumnNamesAcr,levels=c('Anx.','B.S.','C.M.O.P.','C.Im.','Disemb.','E.Im.','E.U.','I.C.','Ins.','S.E.','Syn.'))

# bar version
ggplot(df, aes(x = Values, y = ColumnNamesAcr,fill=Sig)) +
  geom_bar(size=3,stat = "identity")+
  scale_fill_manual(values = colors) +
  labs(y = "", x = "Subscale Correlation",title="Change in Prop. Direction")+theme_minimal(base_size=17)+xlim(c(-.3,.8))+
  theme(legend.position='none')
```

![](Stats_n_Viz_files/figure-gfm/unnamed-chunk-21-3.png)<!-- -->

``` r
# and one just for the legend
df$Significant=df$Sig
ggplot(df, aes(x = Values, y = ColumnNamesAcr,fill=Significant)) +
  geom_bar(size=3,stat = "identity")+
  scale_fill_manual(values = colors) +
  labs(y = "", x = "Subscale Correlation",title="Change in Prop. Direction")+theme_minimal(base_size=17)+xlim(c(-.3,.8))
```

![](Stats_n_Viz_files/figure-gfm/unnamed-chunk-21-4.png)<!-- -->

``` r
# DMN Mag
# Create a dataframe with the values and column names
df <- data.frame(Values = M_corvec, pvals=M_pvec)

# make a color vector by significance
colorvec=rep('NS',11)
colorvec[M_pvec<.05]='Uncr'
colorvec[M_pvec_mc<.05]='FDR'
df$Sig <- factor(colorvec, levels = c('NS', 'Uncr', 'FDR'))
# Specify color scale manually
colors <- c('NS' = '#002642', 'Uncr' = '#EF9500', 'FDR' = '#840032')
df$ColumnNamesAcr<-c(dascNamesAcr)
# bar version
ggplot(df, aes(x = Values, y = ColumnNamesAcr,fill=Sig)) +
  geom_bar(size=3,stat = "identity")+
  scale_fill_manual(values = colors) +
  labs(y = "", x = "Subscale Correlation",title="Change in Prop. Magnitudes")+theme_minimal(base_size=17)+xlim(c(-.3,.8))+
  theme(legend.position='none')
```

    ## Warning: Removed 1 row containing missing values or values outside the scale range
    ## (`geom_bar()`).

![](Stats_n_Viz_files/figure-gfm/unnamed-chunk-21-5.png)<!-- -->

``` r
# pull out BUP~Impaired Control correlation
ggplot(data=DAS_changeDF_people,aes(y=dascscore_impair,x=BUP_Decrease))+
  geom_point(size=4,aes(color = People))+
  labs(x = "Decrease in Bottom-up %", y = "Increase in Impaired Control", color = "People") +
  scale_color_manual(values = generated_colors) +
  theme_minimal(base_size = 16)
```

![](Stats_n_Viz_files/figure-gfm/unnamed-chunk-21-6.png)<!-- -->

``` r
cor.test(DAS_changeDF_people$dascscore_impair,DAS_changeDF_people$BUP_Decrease)
```

    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  DAS_changeDF_people$dascscore_impair and DAS_changeDF_people$BUP_Decrease
    ## t = 4.4524, df = 20, p-value = 0.0002444
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  0.4041507 0.8687335
    ## sample estimates:
    ##       cor 
    ## 0.7055399
