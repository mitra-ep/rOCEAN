#' @title ocean algorithm 
#'
#' @description Estimates TDP in 3 scales for a given region of interest
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
#' @param saveP If true, saves the correlation matrix
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

ocean<-function(mps, om1, om2, p1, p2, gCT, scale=c("pair","row","col"),
                cor=FALSE, nMax=100){
  
  #check the pvalue matrix
  if(cor==T & missing(mps)) stop("Matrix of p-values not provided.")
  if(cor==T) mps<-as.matrix(mps)
  if(cor==T & c(nrow(mps)<1 | ncol(mps)<1) ) stop("Check the p-value matrix.")
  
  #parameters
  grandH=gCT[1]
  concp=gCT[2]
  z=gCT[3]
  alpha=gCT[4]
  m=gCT[5]
  
  if(missing(scale)){
    scale=c("pair","row","col")
  }
  
  #TD over pairs
  if("pair" %in% scale){
   
    if(missin(mps)){
    #calculate vector of p values
    pps<-corPs(om1, om2, p1, p2, type="Vec",pthresh =concp )
    if(length(pps)==0) ptdp=0 }

    #run pairwise algorithm
    if(!missing(mps) & missing(ptdp)){
    pps<-as.vector(mps)
    ptdp<-pairTD(pps, m, gCT)
    }
    
  }else{ ptdp<-NA }
  
  gc()
  
  if(sum(c("row","col") %in% scale)>=1){
    #calculate matrix of p values
   if(missing(mps)) mps<-corPs(om1, om2, p1, p2, type="Mat")}
    
  if("row" %in% scale){
      sCatr<-getCat(mps, gCT, m, scale="row")
      
      #run single step algorithm
      ssr<-singleStep(sCatr)
      
      #td if conclusive result                    
      rtdp=1-(ssr$bound/nrow(sCatr) )
      stepr=1

      #td if inconclusive and run BaB    
      if(ssr$bound!=ssr$heuristic){ 
          bboutr<-runbab(sCatr,ssr$heuristic,ssr$bound,nMax=nMax)
          stepr=bboutr$Step }
      
        } else{ rtdp=NA
                stepr=NA
         }
  

    if("col" %in% scale){
        sCatc<-getCat(mps, gCT, m, scale="col")
        
        #run single step algorithm
        ssc<-singleStep(sCatc)
        
        #td if conclusive result                    
        ctdp=1-(ssc$bound/nrow(sCatc) )  
        stepc=1

        #td if inconclusive and run BaB    
        if(ssc$bound!=ssc$heuristic){ 
            bboutc<-runbab(sCatc,ssc$heuristic,ssc$bound,nMax=nMax)
            stepc=bboutc$Step }
        
          }  else{ ctdp=NA
                   stepc=NA}
  
    #arrange items to return
    #result for row
    if(!is.na(stepr) & stepr>1){
    outr<-c("rHeuristic"=bboutr$sL,"rBound"=bboutr$sU, "nStep"=stepr)
    }else {outr<-c("rTDP"=rtdp, "nStep"=stepr)}
  
    #result for row
    if(!is.na(stepc) & stepc>1){
    outc<-c("cHeuristic"=bboutc$sL,"cBound"=bboutc$sU, "nStep"=stepc)
    }else {outc<-c("cTDP"=ctdp, "nStep"=stepc)}
  
  
    out<-list("Pair"=c("TDP"=ptdp),
              "Row"=outer,
              "Col"=outc)
      
    ###remove NAs from item to return
    out<-lapply(out, function(x) x[!is.na(x)])
    out<-out[lapply(out,length)>0]
    
    return(out)
  }
