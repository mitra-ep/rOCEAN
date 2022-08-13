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
  concp=gCT[2]
  z=gCT[3]
  alpha=gCT[4]
  m=as.numeric(nrow(om1))*as.numeric(nrow(om2))
  
  if(missing(scale)){
    scale=c("pair","row","col")
  }
  #FD over pairs
  if("pair" %in% scale){
   
    #calculate vector of p values
    cat("Calculating P vector. \n")
    pps<-corPs(om1, om2, p1, p2, type="Vec",pthresh =concp )
    if(length(pps)==0) pfd=m

    #run pairwise algorithm
    if(length(pps)>0){
    cat("Calculating FD for pairs. \n")
    pfd<-pairFD(pps, gCT)}
    
  }else{ pfd<-NA }
  
  gc()
  
  if(sum(c("row","col") %in% scale)>=1){
    #calculate matrix of p values
    ps<-corPs(om1, om2, p1, p2, type="Mat")}
    
  if("row" %in% scale){
      sCatr<-getCat(ps, gCT, m, scale="row")
      
      #run single step algorithm
      ssr<-singleStep(sCatr)
      
      #fd if conclusive result                    
      rfd=ssr$bound-1
      
      #fd if inconclusive and run BB    
      if(BB==TRUE) {
        if(ssr$bound-1!=ssr$heuristic){ 
          bboutr<-runbab(sCatr,ssr$heuristic,ssr$bound)
          brfd=bboutr$FD
          }else brfd=NA
         } 
           } else{ rfd=NA
                   if(BB==TRUE) brfd=NA}
  

    if("col" %in% scale){
        sCatc<-getCat(ps, gCT, m, scale="col")
        
        #run single step algorithm
        ssc<-singleStep(sCatc)
        
        #fd if conclusive result                    
        cfd=ssc$bound-1   
          
        #fd if inconclusive and run BB    
        if(BB==TRUE) {
          if(ssc$bound-1!=ssc$heuristic){ 
            bboutc<-runbab(sCatc,ssc$heuristic,ssc$bound)
            bcfd=bboutc$FD
            }else bcfd=NA
          } 
              } else{ cfd=NA
                      if(BB==TRUE) bcfd=NA}
  
    #arrange items to return
    if(BB==TRUE) {
      out<-list("SingleStep"=c("pairfd"=pfd/m,
                                          "SSr"=rfd/length(p1),
                                          "SSc"=cfd/length(p2) ) ,
                           
                           "BandB"=c("b.fdr"=brfd/length(p1),
                                   "b.fdc"=bcfd/length(p2) ) ) }
      
    if(BB==FALSE) {

        out<-c("SingleStep"=c("pairfd"=pfd/m,
                                 "SSr"=rfd/length(p1),
                                 "SSc"=cfd/length(p2) )  ) }
        
    ###remove NAs from item to return
    out<-lapply(out, function(x) x[!is.na(x)])
    out<-out[lapply(out,length)>0]
    
    return(out)
  }
