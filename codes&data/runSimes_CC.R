#data
#original data
load("/path_to_files/data_exp_colon.RData")
load("/path_to_files/data_CN_colon.RData")

expdat<-data.exp[,-c(1,2,3)]
cndat<-data.snp[,-c(1,2,3)]

#remove original files
rm(data.exp,data.snp)
gc()

funclist<-paste0("/path_to_files/", c("corPs.R","simesCT.R") )
a<-lapply(funclist, source)
st<-Sys.time()
CT<-simesCT(expdat, cndat, alpha=0.05)
et<-Sys.time()
as.numeric(difftime(et, st, units="min"))

save(CT, file = "/path_to_files/gCT_colon.RData")

