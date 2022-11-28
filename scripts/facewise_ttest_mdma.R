# extract range of vertices to be covered in this run from VertBin
Lfaces=4589
Rfaces=4595

### MANUALLY RECREATE SUBJ LIST
subjPrefix=rep('sub-MDMA0',17)
subjSuffix=c('01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17')
subjList=paste0(subjPrefix,subjSuffix)

# initialize face-level vectors: keeping 0s in the slots untouched by this run to verify allocation later
T_L=rep(0,Lfaces)
T_R=rep(0,Rfaces)
p_L=rep(0,Lfaces)
p_R=rep(0,Rfaces)
dif_L=rep(0,Lfaces)
dif_R=rep(0,Rfaces)

# create dataframe for output
df<-data.frame(subjList)
colnames(df)<-'SubjID'

# initialize iterable face value column - to avoid storing this data for multiple faces over this upcoming loop
df$BuProp=rep(0,length(subjList))
df$BuProp_mdma=rep(0,length(subjList))

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

# subjvec: no 4
subjvec=c(1,2,3,5,6,7,8,9,10,11,12,13,14,15,16,17)

# for each left hemi face
for (f in 1:Lfaces){
	print(f)
	# load in BU props iteratively
	for (s in subjvec){
		subj=subjList[s]
		#file paths
		LbaselineFP=paste0('/oak/stanford/groups/leanew1/users/apines/OpFlAngDs/mdma/',subj,'/',subj,'_',p1[s],'_rs_BUTD_L.csv')
		Lp2FP=paste0('/oak/stanford/groups/leanew1/users/apines/OpFlAngDs/mdma/',subj,'/',subj,'_',p2[s],'_rs_BUTD_L.csv')
		Lm1FP=paste0('/oak/stanford/groups/leanew1/users/apines/OpFlAngDs/mdma/',subj,'/',subj,'_',m1[s],'_rs_BUTD_L.csv')
		Lm2FP=paste0('/oak/stanford/groups/leanew1/users/apines/OpFlAngDs/mdma/',subj,'/',subj,'_',m2[s],'_rs_BUTD_L.csv')
		# load in dat data
		Lbaseline=read.csv(LbaselineFP)
		Lp2=read.csv(Lp2FP)
		Lm1=read.csv(Lm1FP)
		Lm2=read.csv(Lm2FP)
	    	# extract rest BuProp
		BuProp_p1=Lbaseline[f,1]
		BuProp_p2=Lp2[f,1]
		# averaging: mean buprop across placebo and baseline
		df$BuProp[s]=mean(BuProp_p1,BuProp_p2)
	    	# extract mdma BuProp
		BuProp_m1=Lm1[f,1]
                BuProp_m2=Lm2[f,1]
		# averaging
		df$BuProp_mdma[s]=mean(BuProp_m1,BuProp_m2)
	}
	# remove null row 
	testdf=df[-c(4),]
	# t test em
	ttestres=t.test(testdf$BuProp,testdf$BuProp_mdma,paired=TRUE) 
	# extract t stat
	T_L[f]=ttestres$statistic
	p_L[f]=ttestres$p.value
	dif_L[f]=ttestres$estimate
}

# for each Right hemi face
for (f in 1:Rfaces){
        print(f)
        # load in BU props iteratively
        for (s in subjvec){
                subj=subjList[s]
                #file paths
                RbaselineFP=paste0('/oak/stanford/groups/leanew1/users/apines/OpFlAngDs/mdma/',subj,'/',subj,'_',p1[s],'_rs_BUTD_R.csv')
                Rp2FP=paste0('/oak/stanford/groups/leanew1/users/apines/OpFlAngDs/mdma/',subj,'/',subj,'_',p2[s],'_rs_BUTD_R.csv')
                Rm1FP=paste0('/oak/stanford/groups/leanew1/users/apines/OpFlAngDs/mdma/',subj,'/',subj,'_',m1[s],'_rs_BUTD_R.csv')
                Rm2FP=paste0('/oak/stanford/groups/leanew1/users/apines/OpFlAngDs/mdma/',subj,'/',subj,'_',m2[s],'_rs_BUTD_R.csv')
                # load in dat data
                Rbaseline=read.csv(RbaselineFP)
                Rp2=read.csv(Rp2FP)
                Rm1=read.csv(Rm1FP)
                Rm2=read.csv(Rm2FP)
                # extract rest BuProp
                BuProp_p1=Rbaseline[f,1]
                BuProp_p2=Rp2[f,1]
                # averaging: mean buprop across placebo and baseline
                df$BuProp[s]=mean(BuProp_p1,BuProp_p2)
                # extract mdma BuProp
                BuProp_m1=Rm1[f,1]
                BuProp_m2=Rm2[f,1]
                # averaging
                df$BuProp_mdma[s]=mean(BuProp_m1,BuProp_m2)
        }
	# remove null row 
        testdf=df[-c(4),]
       # t test em
        ttestres=t.test(testdf$BuProp,testdf$BuProp_mdma,paired=TRUE)
        # extract t stat
        T_R[f]=ttestres$statistic
        p_R[f]=ttestres$p.value
        dif_R[f]=ttestres$estimate
}

# saveout t stat and ps - still needs to be merged with results from other hemi for MC correction
write.csv(T_L,paste0('/oak/stanford/groups/leanew1/users/apines/results/p_vs_mdma_L_ts.csv'))
write.csv(p_L,paste0('/oak/stanford/groups/leanew1/users/apines/results/p_vs_mdma_L_ps.csv'))
write.csv(T_R,paste0('/oak/stanford/groups/leanew1/users/apines/results/p_vs_mdma_R_ts.csv'))
write.csv(p_R,paste0('/oak/stanford/groups/leanew1/users/apines/results/p_vs_mdma_R_ps.csv'))
