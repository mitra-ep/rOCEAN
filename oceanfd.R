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


##simple single step
oceanfd<-function(om1, om2, p1, p2, gCT, scale=c("pair","row","col")){
  
  #parameters
  grandH=gCT[1]
  z=gCT[2]
  alpha=gCT[3]
  m=nrow(om1)*nrow(om2)
  
  if(missing(scale)){
    scale=c("pair","row","col")
  }
  #FD over pairs
  if("pair" %in% scale){
   
    #calculate vector of pvalues
    ps<-corPs(om1, om2, p1, p2, type="Vec", pthresh=alpha)
    
    #run pairwise algorithm
    fd<-pairFD(p, gCT, aplpha)
    
    }else{
    #FD over rows or columns
        #calculate matrix of pvalues
        ps<-corPs(om1, om2, p1, p2, type="Mat")
          
        #get corresponding category matrix
        if("row" %in% scale){
        sCat<-getCat(ps, gCT, m, scale="row") }
          
        if("col" %in% scale){
        sCat<-getCat(ps, gCT, m, scale="col") }
        
        #run single step algorithm
        fd<-rcFD(sCat)
                            
    if(fd$Bo-1==fd$heuristic) fd<-fd$heuristic 
        
    if(fd$Bo-1!=fd$heuristic & BB==TRUE) fd<-BB()
  
    if(fd$Bo-1!=fd$heuristic & BB==FALSE) fd<-fd  
    }
  
  ###items to return
  return(fd)
}
