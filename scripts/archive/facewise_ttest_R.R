# extract range of vertices to be covered in this run from VertBin
Lfaces=4589
Rfaces=4595

# load in subj list
subjList=read.delim('~/ispotSubjs.txt',header=F)

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
df$BuProp=rep(0,nrow(subjList))
df$BuProp_gng=rep(0,nrow(subjList))

# calculcate study/subject start/end points here so it is less confusing within the loops
#ispot1start=1
#ispot1end=nrow(subjListI)
#engStart=1+nrow(subjListI)
#engEnd=(engStart)+nrow(subjListE)-1
#ispot2start=engEnd+1
#ispot2end=nrow(subjList)

# for each left hemi face
for (f in 1:Rfaces){
	print(f)
        
       # load in BU props iteratively
       for (s in 1:nrow(subjList)){
           subj=df$SubjID[s]
           ncfFP_L=paste0('/oak/stanford/groups/leanew1/users/apines/OpFlAngDs/ispot/',subj,'/ncf_BUTD_R.csv')
           # load in dat data
           ncf=read.csv(ncfFP_L)
           gngFP_L=paste0('/oak/stanford/groups/leanew1/users/apines/OpFlAngDs/ispot/',subj,'/gng_BUTD_R.csv')
           # load in dat data
           gng=read.csv(gngFP_L)
           # extract rest BuProp length
           df$BuProp[s]=ncf[f,1]       
           # extract carit BuProp length
           df$BuProp_gng[s]=gng[f,1]
	}

	# load in BU props iteratively
#	for (s in ispot1start:ispot1end){
#	    subj=df$SubjID[s]
#	    ncfFP_L=paste0('/oak/stanford/groups/leanew1/users/apines/task_contrasts/ispot/',subj,'/ncf_BUTD_L.csv')
	    # load in dat data
#	    ncf=read.csv(ncfFP_L)
 #   	    gngFP_L=paste0('/oak/stanford/groups/leanew1/users/apines/task_contrasts/ispot/',subj,'/gng_BUTD_L.csv')
            # load in dat data
#            gng=read.csv(gngFP_L)
	    # extract rest BuProp length
#	    df$BuProp[s]=ncf[f,1]	
	    # extract carit BuProp length
#	    df$BuProp_gng[s]=gng[f,1]
#	}
	# load in engage loops
#        for (s in engStart:engEnd){
#            subj=df$SubjID[s]
#            ncfFP_L=paste0('/oak/stanford/groups/leanew1/users/apines/task_contrasts/engage/',subj,'/ncf_BUTD_L.csv')
            # load in dat data
#            ncf=read.csv(ncfFP_L)
#            gngFP_L=paste0('/oak/stanford/groups/leanew1/users/apines/task_contrasts/engage/',subj,'/gng_BUTD_L.csv')
            # load in dat data
#            gng=read.csv(gngFP_L)
#            # extract rest BuProp length
#            df$BuProp[s]=ncf[f,1]       
#            # extract carit BuProp length
#            df$BuProp_gng[s]=gng[f,1]
#        }
	# now ispot 2 loops
#	for (s in ispot2start:ispot2end){
#            subj=df$SubjID[s]
#            ncfFP_L=paste0('/oak/stanford/groups/leanew1/users/apines/task_contrasts/ispot2/',subj,'/ncf_BUTD_L.csv')
#            # load in dat data
#            ncf=read.csv(ncfFP_L)
#            gngFP_L=paste0('/oak/stanford/groups/leanew1/users/apines/task_contrasts/ispot2/',subj,'/gng_BUTD_L.csv')
#            # load in dat data
#            gng=read.csv(gngFP_L)
#            # extract rest BuProp length
#            df$BuProp[s]=ncf[f,1]       
#            # extract carit BuProp length
#            df$BuProp_gng[s]=gng[f,1]
#        }	
	# t test em
	ttestres=t.test(df$BuProp,df$BuProp_gng,paired=TRUE) 
	# extract t stat
	T_R[f]=ttestres$statistic
	p_R[f]=ttestres$p.value
	dif_R[f]=ttestres$estimate
}

# saveout dr2s and ps - still needs to be merged with results from other hemi for MC correction
write.csv(T_R,paste0('/oak/stanford/groups/leanew1/users/apines/results/ncf_vs_gng_R_ts.csv'))
write.csv(p_R,paste0('/oak/stanford/groups/leanew1/users/apines/results/ncf_vs_gng_R_ps.csv'))
