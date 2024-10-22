#' @title pairwise true discoveries proportion
#'
#' @description Calculates the PTD over pair; based on SEA algorithm
#'
#' @param p  Matrix of pairwise associations, can be passed as a vector
#' 
#' @param n  Number of pairs; may not be the same as number of pairs in p if a
#' threshold is used to remove large p-values
#' 
#' @param gCT Parameters of the global closed testing provided as the output of simesCT function 
#'  
#' @return Proportion of true discoveries out of n pairs
#'
#' @author Mitra Ebrahimpoor
#'
#' \email{m.ebrahimpoor@@lumc.nl}
#'
#' @seealso [simesCT()]
#' 
#' @export
#' 
#' 

pairTD<-function(p,
                 n,
                 gCT){
  
  #parameters
  grandH=gCT[1]
  alpha=gCT[4]
  
  #checl p
  if(length(p)==0){
    d=0
  }else{
  #a vector
  p<-as.vector(p)
  #sort pvals
  sp<-sort(p)
  
  #create vector of U values
  uval<-sapply(1:length(sp), function(u) 1-u+sum( (sp*grandH)<=(u*alpha) ))
  d<-max(uval)}
  
  #tdp
  tdp<-d/n

  ###return
  return(tdp)
}
