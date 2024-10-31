#' @title pairwise true discoveries proportion
#'
#' @description Calculates the TDP over pairs; based on SEA algorithm
#'
#' @param mps  Matrix or vector of pairwise associations.
#' 
#' @param n  Number of pairs; may not be the size of p if a threshold is used to remove large p-values.
#' 
#' @param gCT Parameters of the global closed testing, output of simesCT function. 
#'  
#' @return Proportion of true discoveries out of n pairs of features.
#'
#'
#' @seealso \link[rSEA]{SEA}, \link{simesCT}
#' 
#' @export
#' 
#' 

pairTDP<-function(mps, n, gCT){
  
  #parameters
  grandH=gCT[1]
  alpha=gCT[4]
  
  #donotrun if no signal
  if(gCT[3]==0) stop("No discoveries in this data!")
  
  #make sure vector and thresh applied
  mps<-as.vector(mps[mps<alpha])
  
  #check p-matrix
  if(length(mps)==0){
    d=0
  }else{
  #vector
  mps<-as.vector(mps)
  #sort pvals
  sp<-sort(mps)
  
  #create vector of U values
  uval<-sapply(1:length(sp), function(u) 1-u+sum( (sp*grandH)<=(u*alpha) ))
  d<-max(uval, na.rm=TRUE)}
  
  #tdp
  tdp<-d/n

  ###return
  return(tdp)
}
