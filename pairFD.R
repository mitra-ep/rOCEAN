#' @title pairFD
#'
#' @description Calculate number of false discoveries over pairs
#'
#' @param p  vector of p-values
#' 
#' @param gCT Parameters of global closed testing as returned by getH
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


##simple single step
pairFD<-function(p, gCT, aplpha){
  
  #sort pvals
  sp<-sort(tps)
  #create vector of hp/alpha
  tp<-(grandH*sp)/alpha
  d<-max(sapply(1:length(tp), function(u) 1-u+sum(tp<u) ))
  
  #fdp
  fdp<-d/length(p)
  ###items to return
  return(list("heuristic"=Hu, "bound"=Bo))
}
