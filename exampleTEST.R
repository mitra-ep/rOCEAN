#load all functions
funclist<-c("corPs.R","getCat.R","oceanfd.R","pairFD.R",
            "runBaB.R","simesCT.R","SingleStep.R") 
z<-lapply(funclist, source)

#test the functions from OCEAN
set.seed(250)
pmat<-matrix(runif(7000,0.05,1)^4, nrow = 70, byrow =T )

gh<-hommel::hommel(as.vector(pmat))
d<-hommel::discoveries(gh)
m<-nrow(pmat)*ncol(pmat)
z<-sum(pmat<=hommel::concentration(gh))
gcvec<-c(m-d,z,0.05)

#not improved by bb
xs<-sample(1:nrow(pmat), 12)  
ys<-sample(1:ncol(pmat), 35) 
scat<-getCat(pmat[xs,ys], CT, m, scale="col")
singleStep(scat)
ss1<-singleStep(scat)
bbout<-runbab(scat,ss1$heuristic,ss1$bound,nMax = 100)
bbout$FD
bbout$Bd
#bbout$allout


pairFD(as.vector(pmat[xs,ys]), 12*35, CT)

#improved by BB
x1<-c(18, 31, 55, 12, 44, 47, 39, 14, 15, 13, 25, 38)
y1<-c(37, 31, 88, 23, 17, 97, 89, 72, 42, 86, 95, 83, 35,
      9, 24, 69, 76, 21, 93, 27,  5, 75, 46, 30, 20, 100,
      65, 45, 62, 19, 59, 67, 50, 73, 66)
scat<-getCat(pmat[x1,y1], gcvec, m, scale="row")
singleStep(scat)

ss1<-singleStep(scat)
bbout<-runbab(scat,ss1$heuristic,ss1$bound,nMax = 100)
bbout$FD
bbout$Bd
#bbout$allout