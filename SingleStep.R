#' @title Single step algorithm 
#'
#' @description Calculates a heuristic and an upper-bound for the number of FD
#'
#' @param pathlist  \code{subject}.
#'
#' @return A list of two objects, the heuristic and the bound
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
#' @importFrom ff ff
#' 
#' 


##simple single step
singleStep<-function(sCat){
  
  #calculate cumsum over columns for customized mat
  cumcat.p<-apply(sCat, 2, cumsum)
  #get R for heuristic
  getRp<-apply(cumcat.p, 1, function(b) any(b>=1:ncol(sCat)))
  
  #calculate cumsum over columns for customized mat
  cumcat.s<-apply(apply(sCat, 2, sort), 2, cumsum) 
  #get R for bound
  getRs<-apply(cumcat.s, 1, function(b) any(b>=1:ncol(sCat)))  
  
  ##calculate the heuristic for all cases
  #extreme cases
  if(all(getRp) | all(!getRp) ) {
    #All rows rejected
    if(all(getRp)) Hu<-0
    #no rows rejected
    if(all(!getRp)) Hu<-length(getRp)
    
    #other cases 
  }else Hu<-max(which(!getRp)) 
  
  
  ##calculate the upper-bound for all cases
  #extreme cases
  if(all(getRs) | all(!getRs) ) {
    if(all(getRs)) Bo<-1
    if(all(!getRs)) Bo<-length(getRp)+1
    
    #other cases 
  }else Bo<-min(which(getRs)) 
  
  ###items to return
  return(list("heuristic"=Hu, "bound"=Bo))
}
