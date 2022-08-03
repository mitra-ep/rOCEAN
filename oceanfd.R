#' @title ocean algorithm 
#'
#' @description Estimated FD in 3 scales for a given region of interest, if single step is inconclusive,
#' it is possible to run branch and bound algorithm
#'
#' @param om1,om2  The omics datasets where rows are probes and columns are samples.
#'
#' @param p1,p2  Optional argument to select pathways based on indexes (or row-names) of the omics datasets
#' which define the pathways of interest corresponding to each omic dataset
#'
#' @param gCT Parameters of the global closed testing provided as the output of simesCT function 
#' 
#' @param scale Optional, will specify the scale of the quantification, a character string. Possible choices are "pair", "col" and "row".
#' If not provided, all scales are returned.
#' 
#' @return For pair, fd is returned. For row and col, FD is returned if single step is conclusive 
#'
#' @author Mitra Ebrahimpoor
#'
#' \email{m.ebrahimpoor@@lumc.nl}
#'
#' @seealso
#'
#' @references
#'
#' @examples
#'
#' @export
#' 
#' 

oceanfd<-function(om1, om2, p1, p2, gCT, scale=c("pair","row","col"), BB=TRUE){
  
  #parameters
  grandH=gCT[1]
  z=gCT[2]
  alpha=gCT[3]
  m=as.numeric(nrow(om1))*as.numeric(nrow(om2))
  
  if(missing(scale)){
    scale=c("pair","row","col")
  }
  #FD over pairs
  if("pair" %in% scale){
   
    #calculate vector of p values
    cat("Calculating P vector. \n")
    pps<-corPs(om1, om2, p1, p2, type="Vec", pthresh=0.01)
    
    #run pairwise algorithm
    cat("Calculating FD for pairs. \n")
    pfd<-pairFD(pps, gCT)
    
  }else{ pfd<-NA }
  
  if(sum(c("row","col") %in% scale)>=1){
    #calculate matrix of p values
    ps<-corPs(om1, om2, p1, p2, type="Mat")}
    
  if("row" %in% scale){
      sCatr<-getCat(ps, gCT, m, scale="row")
      
      #run single step algorithm
      ssr<-singleStep(sCatr)
      
      #fd if conclusive result                    
      if(ssr$Bo-1==ssr$heuristic) rfd=ssr$heuristic 
  
      #fd if inconclusive and BB run   
     if(ssr$Bo-1!=ssr$heuristic & BB==TRUE) {
        bboutr<-runbab(sCatr)
        rfd=bboutr$FD
      }
      
      #return inconclusive
      if(ssr$Bo-1!=ssr$heuristic & BB==FALSE) rfd=ssr$heuristic 
     } else{ rfd<-NA }
      
      
    if("col" %in% scale){
        sCatc<-getCat(ps, gCT, m, scale="col")
        
        #run single step algorithm
        ssc<-singleStep(sCatc)
        
        #fd if conclusive result                    
       if(ssc$Bo-1==ssc$heuristic) cfd=ssc$Bo-1   
          
        #fd if inconclusive and BB run   

       if(ssc$Bo-1!=ssc$heuristic & BB==TRUE) {
         bboutc<-runbab(sCatc)
         cfd=bboutc$FD
          } 
          
        #return inconclusive
        if(ssc$Bo-1!=ssc$heuristic & BB==FALSE) cfd=ssc$Bo-1  
         } else{ cfd<-NA }
  
  #arrange items to return
    
    ###items to return
    return(list("fd"=pfd,"SSr"=rfd,"SSc"=cfd))
  }
