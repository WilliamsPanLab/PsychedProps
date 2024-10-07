# needed libraries
library(nlme)
# load in info from R
# will need mean FD for DES

subjInfo=readRDS('/oak/stanford/groups/leanew1/users/apines/forVertexwise_DES.rds')
# convert subjinfo to format of vertexwise csvs
subjInfo=subjInfo[,c('Subjects','Task','FD','RemTRs')]
colnames(subjInfo)[1]='Subject'
colnames(subjInfo)[3]='MeanFD'
subjInfo$Task=as.factor(subjInfo$Task)
subjInfo$Subject=as.factor(subjInfo$Subject)

# initialize vertex-level vectors
TaskT_L=rep(0,2562)
TaskT_R=rep(0,2562)
Taskp_L=rep(0,2562)
Taskp_R=rep(0,2562)

# for each vertex
for (v in 1:2562){
	print(v)
	# if file exists: left
	if (file.exists(paste0('/scratch/users/apines/taskVerts/v',v,'_DES_L.csv'))){
		# load in data for this vertex
		dataV=read.csv(paste0('/scratch/users/apines/taskVerts/v',v,'_DES_L.csv'))
	
		# one value per scan
		for (s in 1:length(unique(dataV$Subject)){
		     curSubj=unique(dataV$Subject)[s]
		     # get mean value per 
		}
		# remove every other opfl measurement so no single TR is used twice in observations
		
		dataV <- dataV[seq(2, nrow(dataV), by = 2), ]
		# convert task
		dataV$Task[dataV$Task=='rs1']='rs'
		# combine with subjinfo
		combinedData=merge(dataV,subjInfo,by=c('Subject','Task'))
		# try averaging every value for each task for each subject:

		# temp: test if there are any rows in combined data not represented in subjInfo (those that passed QC)
		#combinedData$key <- paste(combinedData$Subject, combinedData$Task, combinedData$Session, sep = "_")
		#subjInfo$key <- paste(subjInfo$Subject, subjInfo$Task, subjInfo$Session, sep = "_")
		# Find keys in combinedData that are not in subjInfo
		#mismatchedKeys <- setdiff(combinedData$key, subjInfo$key)
		# Subset the mismatched rows
		#mismatchedRows <- combinedData[combinedData$key %in% mismatchedKeys, ]
		# set rs1 and rs2 to equivalent
		combinedData$Task[combinedData$Task=='rs2']='rs'
		combinedData$Task<-as.factor(combinedData$Task)
		# fit model
		model <- lme(Value ~ MeanFD + Task + RemTRs, random = ~ 1 | Subject, data = combinedData)
		modeltable=summary(model)$tTable
		# print out stats
		TaskT_L[v]=modeltable['Taskwm','t-value']
		Taskp_L[v]=modeltable['Taskwm','p-value']
	}
	# if right file exists
	if (file.exists(paste0('/scratch/users/apines/taskVerts/v',v,'_R.csv'))){
                # load in data for this vertex
                dataV=read.csv(paste0('/scratch/users/apines/taskVerts/v',v,'_R.csv'))
		# remove every other opfl measurement so no single TR is used twice in observations
		dataV <- dataV[seq(2, nrow(dataV), by = 2), ]
                # convert task
                dataV$Task[dataV$Task=='rs1']='rs'
                # combine with subjinfo
                combinedData=merge(dataV,subjInfo,by=c('Subject','Task','Session'))
                combinedData$Task[combinedData$Task=='rs2']='rs'
		combinedData$Task<-as.factor(combinedData$Task)
                # fit model
                model <- lme(Value ~ MeanFD + Drug+Task+Drug*Task + RemTRs, random = ~ 1 | Subject, data = combinedData)
                modeltable=summary(model)$tTable
                # print out stats
                DrugT_R[v]=modeltable['Drug1','t-value']
                Drugp_R[v]=modeltable['Drug1','p-value']
                TaskT_R[v]=modeltable['Taskwm','t-value']
                Taskp_R[v]=modeltable['Taskwm','p-value']
                DrugTaskT_R[v]=modeltable['Drug1:Taskwm','t-value']
                DrugTaskp_R[v]=modeltable['Drug1:Taskwm','p-value']
	}
# end for each vertex
}
# pull out p's where T!=0 : left 
validVerts_L = which(DrugT_L != 0)
valid_Drug_Ps_L=Drugp_L[validVerts_L]
valid_Task_Ps_L=Taskp_L[validVerts_L]
valid_DrugTask_Ps_L=DrugTaskp_L[validVerts_L]
# pull out p's where T!=0 : right
validVerts_R = which(DrugT_R != 0)
valid_Drug_Ps_R=Drugp_R[validVerts_R]
valid_Task_Ps_R=Taskp_R[validVerts_R]
valid_DrugTask_Ps_R=DrugTaskp_R[validVerts_R]
# mc correct
valid_Drug_Ps=c(valid_Drug_Ps_L,valid_Drug_Ps_R)
valid_Task_Ps=c(valid_Task_Ps_L,valid_Task_Ps_R)
valid_DrugTask_Ps=c(valid_DrugTask_Ps_L,valid_DrugTask_Ps_R)
mc_Drug_Ps=p.adjust(valid_Drug_Ps,method='fdr')
mc_Task_Ps=p.adjust(valid_Task_Ps,method='fdr')
mc_DrugTask_Ps=p.adjust(valid_DrugTask_Ps,method='fdr')
# parse back to left and right
mc_Drug_Ps_L=mc_Drug_Ps[1:length(validVerts_L)]
mc_Drug_Ps_R=mc_Drug_Ps[(length(validVerts_L)+1):(length(validVerts_L)+length(validVerts_R))]
mc_Task_Ps_L=mc_Task_Ps[1:length(validVerts_L)]
mc_Task_Ps_R=mc_Task_Ps[(length(validVerts_L)+1):(length(validVerts_L)+length(validVerts_R))]
mc_DrugTask_Ps_L=mc_DrugTask_Ps[1:length(validVerts_L)]
mc_DrugTask_Ps_R=mc_DrugTask_Ps[(length(validVerts_L)+1):(length(validVerts_L)+length(validVerts_R))]
# apply to t's
valid_Drug_Ts_L=DrugT_L[validVerts_L]
valid_Task_Ts_L=TaskT_L[validVerts_L]
valid_DrugTask_Ts_L=DrugTaskT_L[validVerts_L]
# pull out p's where T!=0 : right
valid_Drug_Ts_R=DrugT_R[validVerts_R]
valid_Task_Ts_R=TaskT_R[validVerts_R]
valid_DrugTask_Ts_R=DrugTaskT_R[validVerts_R]
# use FDR < .05 to thresh: inputting 999 for now to track through vis
valid_Drug_Ts_L[mc_Drug_Ps_L>.05]=999
valid_Drug_Ts_R[mc_Drug_Ps_R>.05]=999
valid_Task_Ts_L[mc_Task_Ps_L>.05]=999
valid_Task_Ts_R[mc_Task_Ps_R>.05]=999
valid_DrugTask_Ts_L[mc_DrugTask_Ps_L>.05]=999
valid_DrugTask_Ts_R[mc_DrugTask_Ps_R>.05]=999
# saveout corrected t's for vis
write.csv(data.frame(valid_Drug_Ts_L), file = "/scratch/users/apines/taskVerts/valid_Drug_Ts_L.csv",row.names = FALSE)
write.csv(data.frame(valid_Drug_Ts_R), file = "/scratch/users/apines/taskVerts/valid_Drug_Ts_R.csv",row.names = FALSE)
write.csv(data.frame(valid_Task_Ts_L), file = "/scratch/users/apines/taskVerts/valid_Task_Ts_L.csv",row.names = FALSE)
write.csv(data.frame(valid_Task_Ts_R), file = "/scratch/users/apines/taskVerts/valid_Task_Ts_R.csv",row.names = FALSE)
write.csv(data.frame(valid_DrugTask_Ts_L), file = "/scratch/users/apines/taskVerts/valid_DrugTask_Ts_L.csv",row.names = FALSE)
write.csv(data.frame(valid_DrugTask_Ts_R), file = "/scratch/users/apines/taskVerts/valid_DrugTask_Ts_R.csv",row.names = FALSE)

