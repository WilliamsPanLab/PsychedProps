Lfaces=1800
Rfaces=1898
# initialize full-length DMN L vector
DMNL=rep(0,Lfaces);
DMNLp=rep(0,Lfaces);
# initialize full-length DMN R vector
DMNR=rep(0,Rfaces);
DMNRp=rep(0,Rfaces);
# read in each of 6 /scratch/users/apines/DMN/DMN_Prop_Tstats_L_Drug_',as.character(iteration),'.csv'
for (i in 1:6){
	start_iteration <- (i - 1) * 300 + 1
	end_iteration <- i * 300
	fp=paste0('/scratch/users/apines/DMN/DMN_Prop_Tstats_L_Drug_',as.character(i),'.csv')
	data=read.csv(fp)
	DMNL[start_iteration:end_iteration]=data[start_iteration:end_iteration,2]
	# p
	fp=paste0('/scratch/users/apines/DMN/DMN_Prop_ps_L_Drug_',as.character(i),'.csv')
        data=read.csv(fp)
        DMNLp[start_iteration:end_iteration]=data[start_iteration:end_iteration,2]
}
# read in each of 8 /scratch/users/apines/DMN/DMN_Prop_Tstats_R_Drug_',as.character(iteration),'.csv'
for (i in 1:8){
        start_iteration <- (i - 1) * 271 + 1
        end_iteration <- i * 271
        # this is dumb, but its a catch for the last right-hemisphere face
	if(i==8){
	        start_iteration=1898
	        end_iteration=1898
	}
	fp=paste0('/oak/stanford/groups/leanew1/users/apines/results/DMNProps_Tstats_R_Drug_',as.character(i),'.csv')
        data=read.csv(fp)
        DMNR[start_iteration:end_iteration]=data[start_iteration:end_iteration,2]
	# p
        fp=paste0('/oak/stanford/groups/leanew1/users/apines/results/DMNProps_ps_R_Drug_',as.character(i),'.csv')
        data=read.csv(fp)
        DMNRp[start_iteration:end_iteration]=data[start_iteration:end_iteration,2]
}
# Combine vectors
tsBoth=c(DMNL,DMNR)
psBoth=c(DMNLp,DMNRp)
# FDR correct it all
fdrp=p.adjust(psBoth,method='fdr')
# threshold ts
tsBoth[fdrp>.05]=0
# sep. out left and right again
tL=tsBoth[1:1800]
tR=tsBoth[1801:length(tsBoth)]
# reconstruct into original face-wise mask
OGf_L=read.csv('~/MasterMask_L.csv',header=F)
OGf_R=read.csv('~/MasterMask_R.csv',header=F)
Outvec_L=rep(0,5120)
Outvec_R=rep(0,5120)
Outvec_L[OGf_L==1]=tL
Outvec_R[OGf_R==1]=tR
# save out FDR-corrected reconstructed t's
write.table(Outvec_L,'/scratch/users/apines/gp/PropFeats/Drug_Ts_L_FDRed.csv',col.names=F,row.names=F,quote=F)
write.table(Outvec_R,'/scratch/users/apines/gp/PropFeats/Drug_Ts_R_FDRed.csv',col.names=F,row.names=F,quote=F)
