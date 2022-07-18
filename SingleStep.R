#' @title Single step algorithm 
#'
#' @description Calculates a heuristic and an upper-bound for the number of FD based on the algorithm
#' introduced in paper
#'
#' @param sCat  Category matrix, output of getCat function
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
#' 

singleStep<-function(sCat, B){
  
  ##set B if missing
  if(missing(B)) B=rep(1,nrow(sCat))
  
  
  #B is full-size
  if(length(B)==nrow(sCat)) {
    #one set selected
    if(sum(B)==1)  type="single"
    #all sets mentioned in B
    if(sum(B)>1) type="full"
  }
  
  #B is not full-size
  if(length(B)<nrow(sCat)){
    #one set selected
    if(sum(B)==0 & nrow(sCat)-length(B)==1){
      type="single"}else
        type="partial"
  }
  
  
  #only one single row selected
  if(type=="single"){
    
    submat<-sCat[which(B==1),]
    #get R for heuristic
    getRp<-any(submat>=1:ncol(sCat))
    Hu<-ifelse(getRp,1,0)  
    
    #get R for bound
    getRs<-any(submat>=1:ncol(sCat))
    Bo<-ifelse(getRs,1,0)     }  
  
  #a subset of rows selected and some not changed
  if(type=="partial"){
    #subscripts where no action required
    Bfix<-c(B, rep(2,nrow(sCat)-length(B) ) )
    
    #create a new matrix by removing rows based on B
    submat<-rbind(sCat[which(B==1),], sCat[which(Bfix==2),])
    
    if(sum(Bfix==2)>1) submats<-rbind(sCat[which(B==1),], apply(sCat[which(Bfix==2),], 2, sort))
    if(sum(Bfix==2)==1) submats<-rbind(sCat[which(B==1),], sCat[which(Bfix==2),])
    
    #calculate cumsum over columns for customized mat
    cumcat.p<-apply(submat, 2, cumsum)
    #get R for heuristic
    getRp<-apply(cumcat.p, 1, function(b) any(b>=1:ncol(sCat)))
    
    #calculate cumsum over columns for customized mat
    cumcat.s<-apply(apply(submat, 2, sort), 2, cumsum) 
    #get R for bound
    getRs<-apply(cumcat.s, 1, function(b) any(b>=1:ncol(sCat)))
    
  }
  
  #all rows are indexed in B
  if(type=="full"){
    
    #create a new matrix by removing rows based on B
    submat<-sCat[which(B==1),]

    #calculate cumsum over columns for customized mat
    cumcat.p<-apply(submat, 2, cumsum)
    #get R for heuristic
    getRp<-apply(cumcat.p, 1, function(b) any(b>=1:ncol(sCat)))
    
    #calculate cumsum over columns for customized mat
    cumcat.s<-apply(apply(submat, 2, sort), 2, cumsum) 
    #get R for bound
    getRs<-apply(cumcat.s, 1, function(b) any(b>=1:ncol(sCat)))  
  }
  
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
