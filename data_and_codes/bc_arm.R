#data
#original data
load("/path_to_files/data_exp_breast.RData")
load("/ath_to_files/data_CN_breast.RData")

#remove names 
expdat<-data.exp[,-c(1,2,3)]
rownames(expdat)<-data.exp[,1]
cndat<-data.snp[,-c(1,2,3)]
rownames(cndat)<-data.snp[,1]

#remove original files
rm(data.exp,data.snp)
gc()

#load functions from package
library(OCEAN)

#load pathlist
load("/path_to_files/allarm.RData")
load("/path_to_files/pathHall.RData")
load("/path_to_filest/pathArm.RData")

#run ocean
expids<-which(rownames(expdat) %in% hallpath[[allarm$pexp[r,]]]) 
cnids<-which(rownames(cndat) %in% charmlist[[allarm$pcn[r,]]])

#load gCT results
load("/path_to_files/gCT.RData")
st<-Sys.time()
out<-ocean(om1=expdat[expids,], om2=cndat[cnids,], 
             gCT, scale=c("pair","row","col"))
et<-Sys.time()
as.numeric(difftime(et, st, units="min"))
#get the pathlist
save(out, file = "/path_to_files/qsub_bcarm/obcr,.RData")

