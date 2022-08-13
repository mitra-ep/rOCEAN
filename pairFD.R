#' @title pairFD
#'
#' @description Calculate number of false discoveries over pairs
#'
#' @param p  vector of p-values in the region of interest
#' 
#' @param gCT Parameters of the global closed testing provided as the output of simesCT function 
#'  
#' @return Number of false discoveries among pairs
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

pairFD<-function(p, gCT){
  
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
  
  #fdp
  fd<-length(p)-d

  ###return
  return(fd)
}
