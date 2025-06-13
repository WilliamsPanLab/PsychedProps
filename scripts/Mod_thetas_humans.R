################ libraries
library(nlme)
library(bpngreg)
################ load in MDMA
rs1bv=read.csv('/oak/stanford/groups/leanew1/users/apines/data/rs1_thetas_L_bv.csv',header=F)
rs2bv=read.csv('/oak/stanford/groups/leanew1/users/apines/data/rs2_thetas_L_bv.csv',header=F)
emobv=read.csv('/oak/stanford/groups/leanew1/users/apines/data/emotion_thetas_L_bv.csv',header=F)
gamblingbv=read.csv('/oak/stanford/groups/leanew1/users/apines/data/gambling_thetas_L_bv.csv',header=F)
wmbv=read.csv('/oak/stanford/groups/leanew1/users/apines/data/wm_thetas_L_bv.csv',header=F)
rs1bv$Task='rs'
rs2bv$Task='rs2'
emobv$Task='emotion'
gamblingbv$Task='gambling'
wmbv$Task='wm'
# remove subject 4
rs1bv=rs1bv[-c(4),]
rs2bv=rs2bv[-c(4),]

# load in LSD

# load in psil

# model each face

# save out stats
