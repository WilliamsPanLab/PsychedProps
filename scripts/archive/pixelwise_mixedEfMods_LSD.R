# needed libraries
library(lme4)
library(bigmemory)
library(grDevices)
library(lmerTest)

# initialize vertex-level vectors
DrugTs=array(0,dim=c(67,70))
Drugps=array(0,dim=c(67,70))

# loop over each mouse, load in pixelwise
mList <- c('m2000', 'm7507', 'm7520', 'm7522', 'm7589', 'm7594')
# set dimensions
n1 <- 67
n2 <- 70
n3 <- 892
n4 <- length(mList)
n5 <- 6  # Number of timepoints
faceMatrix5D <- array(0, dim = c(n1, n2, n3, n4, n5))
# populate a 5D array: 67 x 70 x 892 x 6 (6 mice), 6 timepoints
# Loop over each mouse and each timepoint to load the data
for (i in 1:length(mList)) {
  subj <- mList[i]
  for (sesh in 1:n5) {
    # account for rogue missing session mouse
    if (subj=='m7507' && sesh==6){
	faceMatrix5D[,,,i,sesh] = array(NA,dim=c(n1,n2,n3))
	} else {
    # Construct the file path
    file_path <- paste0('/scratch/users/apines/gp/PropFeats/', subj, '_', sesh, '_faceMatrix_gro_pixelwise.csv')
    
    # Read the CSV file
    faceMatrix_reshaped <- read.big.matrix(file_path,type='double')
    
    # Convert data frame to matrix
    faceMatrix_reshaped_matrix <- as.matrix(faceMatrix_reshaped)
    
    # Reshape the matrix back to the original dimensions
    faceMatrix <- array(as.numeric(faceMatrix_reshaped_matrix), dim = c(n1, n2, n3))
    
    # Populate the 5D array
    faceMatrix5D[,,,i,sesh] <- faceMatrix
      }
   }
}

# test dummy matrix
# was done to ensure reshaping matrices is working as intended
#dummyMatrixfp='/scratch/users/apines/gp/PropFeats/dummyMatrix.csv'
#dummMat=read.big.matrix(dummyMatrixfp,type='double')
# reshape
#n1 <- 67
#n2 <- 70
#n3 <- 892
#dummMat<-as.matrix(dummMat)
#dummyMatrix <- array(as.numeric(dummMat), dim = c(n1, n2, n3))
# reconstruction works on dummy matrix!

# for each pixel x
for (x in 1:67){
	print(x)
	for (y in 1:70){
		# get all sober values from mouse 1
		# every other value to not utilize observations derived in part from a single frame twice
		m1SobVals=faceMatrix5D[x,y,seq(1, dim(faceMatrix5D)[3], by = 2),1,1]
		m1SobData=data.frame(Value=m1SobVals, MouseID = rep(mList[1],length(m1SobVals)), Drug = rep(0,length(m1SobVals)))
		m2SobVals=faceMatrix5D[x,y,seq(1, dim(faceMatrix5D)[3], by = 2),2,1]
		m2SobData=data.frame(Value=m2SobVals, MouseID = rep(mList[2],length(m2SobVals)), Drug = rep(0,length(m2SobVals)))
		m3SobVals=faceMatrix5D[x,y,seq(1, dim(faceMatrix5D)[3], by = 2),3,1]
		m3SobData=data.frame(Value=m3SobVals, MouseID = rep(mList[3],length(m3SobVals)), Drug = rep(0,length(m3SobVals)))
		m4SobVals=faceMatrix5D[x,y,seq(1, dim(faceMatrix5D)[3], by = 2),4,1]
		m4SobData=data.frame(Value=m4SobVals, MouseID = rep(mList[4],length(m4SobVals)), Drug = rep(0,length(m4SobVals)))
		m5SobVals=faceMatrix5D[x,y,seq(1, dim(faceMatrix5D)[3], by = 2),5,1]
		m5SobData=data.frame(Value=m5SobVals, MouseID = rep(mList[5],length(m5SobVals)), Drug = rep(0,length(m5SobVals)))
		m6SobVals=faceMatrix5D[x,y,seq(1, dim(faceMatrix5D)[3], by = 2),6,1] 	 
		m6SobData=data.frame(Value=m6SobVals, MouseID = rep(mList[6],length(m6SobVals)), Drug = rep(0,length(m6SobVals)))
		# pull out all drug values
		m1DrugVals=array(faceMatrix5D[x,y,seq(1, dim(faceMatrix5D)[3], by = 2),1,2:6])
                m1DrugData=data.frame(Value=m1DrugVals, MouseID = rep(mList[1],length(m1DrugVals)), Drug = rep(1,length(m1DrugVals)))
                m2DrugVals=array(faceMatrix5D[x,y,seq(1, dim(faceMatrix5D)[3], by = 2),2,2:6])
                m2DrugData=data.frame(Value=m2DrugVals, MouseID = rep(mList[2],length(m2DrugVals)), Drug = rep(1,length(m2DrugVals)))
                m3DrugVals=array(faceMatrix5D[x,y,seq(1, dim(faceMatrix5D)[3], by = 2),3,2:6])
                m3DrugData=data.frame(Value=m3DrugVals, MouseID = rep(mList[3],length(m3DrugVals)), Drug = rep(1,length(m3DrugVals)))
                m4DrugVals=array(faceMatrix5D[x,y,seq(1, dim(faceMatrix5D)[3], by = 2),4,2:6])
                m4DrugData=data.frame(Value=m4DrugVals, MouseID = rep(mList[5],length(m4DrugVals)), Drug = rep(1,length(m4DrugVals)))
                m5DrugVals=array(faceMatrix5D[x,y,seq(1, dim(faceMatrix5D)[3], by = 2),5,2:6])
                m5DrugData=data.frame(Value=m5DrugVals, MouseID = rep(mList[6],length(m5DrugVals)), Drug = rep(1,length(m5DrugVals)))
                m6DrugVals=array(faceMatrix5D[x,y,seq(1, dim(faceMatrix5D)[3], by = 2),6,2:6])
                m6DrugData=data.frame(Value=m6DrugVals, MouseID = rep(mList[6],length(m6DrugVals)), Drug = rep(1,length(m6DrugVals)))
		# combine data
		DataCombined=rbind(m1SobData,m2SobData,m3SobData,m4SobData,m5SobData,m6SobData,m1DrugData,m2DrugData,m3DrugData,m4DrugData,m5DrugData,m6DrugData)
		# omit NAs (scan 6 mouse 7507)
		DataCombined=na.omit(DataCombined)
		DataCombined$MouseID<-as.factor(DataCombined$MouseID)
		# fit model
		model <- lmer(Value ~ Drug + (1 | MouseID), data = DataCombined)
		modeltable=summary(model)$coefficients
		# print out stats
		DrugTs[x,y]=modeltable['Drug','t value']
		Drugps[x,y]=modeltable['Drug','Pr(>|t|)']
		# end for each y pixel
	}
	# end for each x pixel
}

# Apply mask to return valid pixels	
Mask=read.csv('~/MouseMaskBool_67x70.csv',header=F)	
# pull out p's where T!=0 : left 
validPixels = which(Mask != 0)
valid_Drug_Ps=Drugps[validPixels]
valid_Drug_Ts=DrugTs[validPixels]

# mc correct
mc_Drug_Ps=p.adjust(valid_Drug_Ps,method='fdr')
# use FDR < .05 to thresh: inputting 999 for now to track through vis
valid_Drug_Ts[mc_Drug_Ps>.05]=999
# reconstruct in original mask space
outputCsv=array(0,dim=c(67,70))
outputCsv[validPixels]=valid_Drug_Ts
# saveout corrected t's for vis
write.csv(outputCsv, file = "/scratch/users/apines/taskVerts/valid_Drug_Ts_pixels.csv",row.names = FALSE)

