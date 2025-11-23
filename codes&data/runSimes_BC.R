#data
#original data
load("/path_to_files/data_exp_breast.RData")
load("/path_to_files/data_CN_breast.RData")

expdat<-data.exp[,-c(1,2,3)]
cndat<-data.snp[,-c(1,2,3)]

#remove original files
rm(data.exp,data.snp)
gc()

#load functions from package
library(OCEAN)

st<-Sys.time()
CT<-simesCT(expdat, cndat, alpha=0.05)
et<-Sys.time()
as.numeric(difftime(et, st, units="min"))

save(CT, file = "/path_to_files/CT_colon.RData")

