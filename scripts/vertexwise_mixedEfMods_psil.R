# needed libraries
library(lme4)
library(lmerTest)
# load in info from R
subjInfo=readRDS('/oak/stanford/groups/leanew1/users/apines/forVertexwise_psil.rds')
# convert subjinfo to format of vertexwise csvs
subjInfo=subjInfo[,c('Subjects','Task','Session','Drug','FD','RemTRs')]
colnames(subjInfo)[1]='Subject'
subjInfo$Task=as.factor(subjInfo$Task)
subjInfo$Subject=as.factor(subjInfo$Subject)

# initialize vertex-level vectors
DrugT_L=rep(0,2562)
DrugT_R=rep(0,2562)
Drugp_L=rep(0,2562)
Drugp_R=rep(0,2562)
TaskT_L=rep(0,2562)
TaskT_R=rep(0,2562)
Taskp_L=rep(0,2562)
Taskp_R=rep(0,2562)
DrugTaskT_L=rep(0,2562)
DrugTaskp_L=rep(0,2562)
DrugTaskT_R=rep(0,2562)
DrugTaskp_R=rep(0,2562)
# for each vertex
for (v in 1:2562){
	print(v)
	# if file exists: left
	if (file.exists(paste0('/scratch/users/apines/taskVerts/v',v,'_psil_L.csv'))){
		# load in data for this vertex
		dataV=read.csv(paste0('/scratch/users/apines/taskVerts/v',v,'_psil_L.csv'))
		# remove every other opfl measurement so no single TR is used twice in observations
		# not applicable for single-value-per-scan measures
		#dataV <- dataV[seq(2, nrow(dataV), by = 2), ]
		# combine with subjinfo
		combinedData=merge(dataV,subjInfo,by=c('Subject','Task','Session'))
		# temp: test if there are any rows in combined data not represented in subjInfo (those that passed QC)
		#combinedData$key <- paste(combinedData$Subject, combinedData$Task, combinedData$Session, sep = "_")
		#subjInfo$key <- paste(subjInfo$Subject, subjInfo$Task, subjInfo$Session, sep = "_")
		# Find keys in combinedData that are not in subjInfo
		#mismatchedKeys <- setdiff(combinedData$key, subjInfo$key)
		# Subset the mismatched rows
		#mismatchedRows <- combinedData[combinedData$key %in% mismatchedKeys, ]
		combinedData <- within(combinedData, Drug <- relevel(Drug, ref = 2))
		# fit model
		model <- lmer(Value ~ FD + Drug + RemTRs + (1 | Subject), data = combinedData)
		modeltable=summary(model)$coefficients
		# print out stats
		DrugT_L[v]=modeltable['DrugPsilo','t value']
		Drugp_L[v]=modeltable['DrugPsilo','Pr(>|t|)']
	}
	# if right file exists
	if (file.exists(paste0('/scratch/users/apines/taskVerts/v',v,'_psil_R.csv'))){
                # load in data for this vertex
                dataV=read.csv(paste0('/scratch/users/apines/taskVerts/v',v,'_psil_R.csv'))
                # combine with subjinfo
                combinedData=merge(dataV,subjInfo,by=c('Subject','Task','Session'))
		combinedData <- within(combinedData, Drug <- relevel(Drug, ref = 2))
                # fit model
                model <- lmer(Value ~ FD + Drug + RemTRs + (1| Subject), data = combinedData)
                modeltable=summary(model)$coefficients
                # print out stats
                DrugT_R[v]=modeltable['DrugPsilo','t value']
                Drugp_R[v]=modeltable['DrugPsilo','Pr(>|t|)']
	}
# end for each vertex
}
# pull out p's where T!=0 : left 
validVerts_L = which(DrugT_L != 0)
valid_Drug_Ps_L=Drugp_L[validVerts_L]
# pull out p's where T!=0 : right
validVerts_R = which(DrugT_R != 0)
valid_Drug_Ps_R=Drugp_R[validVerts_R]
# mc correct
valid_Drug_Ps=c(valid_Drug_Ps_L,valid_Drug_Ps_R)
mc_Drug_Ps=p.adjust(valid_Drug_Ps,method='fdr')
# parse back to left and right
mc_Drug_Ps_L=mc_Drug_Ps[1:length(validVerts_L)]
mc_Drug_Ps_R=mc_Drug_Ps[(length(validVerts_L)+1):(length(validVerts_L)+length(validVerts_R))]
# apply to t's
valid_Drug_Ts_L=DrugT_L[validVerts_L]
# pull out p's where T!=0 : right
valid_Drug_Ts_R=DrugT_R[validVerts_R]
# use FDR < .05 to thresh: inputting 999 for now to track through vis
valid_Drug_Ts_L[mc_Drug_Ps_L>.5]=999
valid_Drug_Ts_R[mc_Drug_Ps_R>.5]=999
# saveout corrected t's for vis
write.csv(data.frame(valid_Drug_Ts_L), file = "/scratch/users/apines/taskVerts/valid_Drug_Ts_L_psil.csv",row.names = FALSE)
write.csv(data.frame(valid_Drug_Ts_R), file = "/scratch/users/apines/taskVerts/valid_Drug_Ts_R_psil.csv",row.names = FALSE)

