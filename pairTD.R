#' @title pairTD
#'
#' @description Calculate PTD over pairs
#'
#' @param p  vector of p-values in the region of interest
#' 
#' @param n  number of pairs
#' 
#' @param gCT Parameters of the global closed testing provided as the output of simesCT function 
#'  
#' @return Proportion of true discoveries among pairs
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

pairTD<-function(p, n, gCT){
  
  #parameters
  grandH=gCT[1]
  alpha=gCT[4]

  #sort pvals
  sp<-sort(p)
  
  #create vector of hp/alpha
  tp<-(grandH*sp)/alpha
  u<-1:length(tp)
  uval<-1-u+sum(tp<u)
  d<-max(uval)
  
  #tdp
  tdp<-d/n

  ###return
  return(tdp)
}
