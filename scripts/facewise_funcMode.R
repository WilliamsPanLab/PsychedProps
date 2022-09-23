
# load in subj list
subjList=read.delim('~/ispotSubjs.txt',header=F)

# initialize face-level vectors: keeping 0s in the slots untouched by this run to verify allocation later
T_L=rep(0,Lfaces)
T_R=rep(0,Rfaces)
p_L=rep(0,Lfaces)
p_R=rep(0,Rfaces)

# create dataframe for output
df<-data.frame(subjList)
colnames(df)<-'SubjID'

# initialize iterable face value column - to avoid storing this data for multiple faces over this upcoming loop
df$BuProp_Shift=rep(0,nrow(subjList))
df$Activation_Shift=rep(0,nrow(subjList))

# load in BU props and oxyhemoglobin signal iteratively
for (s in 1:nrow(subjList)){
           subj=df$SubjID[s]
	# LEFT - BU PROP
           ncfFP_L=paste0('/oak/stanford/groups/leanew1/users/apines/OpFlAngDs/ispot/',subj,'/ncf_BUTD_L.csv')
           # load in dat data
           ncf=read.csv(ncfFP_L)
           gngFP_L=paste0('/oak/stanford/groups/leanew1/users/apines/OpFlAngDs/ispot/',subj,'/gng_BUTD_L.csv')
           # load in dat data
           gng=read.csv(gngFP_L)
           # extract rest BuProp length
           differenceBU=ncf[,1]-gng[,1]       
           # measure evoked difference
	   df$BuProp_Shift[s]=sum(abs(differenceBU))
	# RIGHT - BU PROP
           ncfFP_R=paste0('/oak/stanford/groups/leanew1/users/apines/OpFlAngDs/ispot/',subj,'/ncf_BUTD_R.csv')
           # load in dat data
           ncf=read.csv(ncfFP_R)
           gngFP_R=paste0('/oak/stanford/groups/leanew1/users/apines/OpFlAngDs/ispot/',subj,'/gng_BUTD_R.csv')
           # load in dat data
           gng=read.csv(gngFP_R)
           # extract rest BuProp length
           differenceBU=ncf[,1]-gng[,1]
           # measure evoked difference
           df$BuProp_Shift[s]=df$BuProp_Shift[s]+sum(abs(differenceBU))
	# LEFT - OXYHEME DIFF
           difFP_L=paste0('/oak/stanford/groups/leanew1/users/apines/OpFlAngDs/ispot/',subj,'/GNG-NCF_L.csv')
           # load in dat data
           difL=read.csv(difFP_L)
           # measure evoked difference
           df$Activation_Shift[s]=sum(abs(difL))	
	# RIGHT - OXYHEME DIFF
           difFP_R=paste0('/oak/stanford/groups/leanew1/users/apines/OpFlAngDs/ispot/',subj,'/GNG-NCF_R.csv')
           # load in dat data
           difR=read.csv(difFP_R)
           # measure evoked difference
           df$Activation_Shift[s]=sum(abs(difR))
	print(s)
}

# saveout dataframe
saveRDS(df,'~/PropShift_ActShift.rds')
