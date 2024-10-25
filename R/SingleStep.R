#' @title Single step algorithm 
#'
#' @description Calculates heuristic and upper-bound for the number of true discoveries based on
#' the Algorithm 1 introduced in paper.
#'
#' @param sCat p-categories matrix, output of getCat function.
#' 
#' @param B Optional, to identify rows to be fixed (1) or removed (0) while splitting the search space.
#'  
#' @return A list of two objects, the heuristic and the lower bound for true number of discoveries
#'
#' @author Mitra Ebrahimpoor
#'
#' \email{m.ebrahimpoor@@lumc.nl}
#'
#' @seealso \link{getCat}
#'
#' @export
#' 
#' 

singleStep<-function(sCat, B){
  
  ##search for the smallest j where v_jk>=k
  findj<-function(mat){
    r<- nrow(mat)
    l<- ncol(mat)
    
    #if no mat_ij>=k
    res=r+1
    
    #initiate while loop
    j<-r
    k<-1
    while (k<=l && j>=1) {
      if(mat[j,k]>=k) {
        res=j
        j<-j-1
      } else {
        k<-k+1
      }
    }
    return(res)
  }
  
  ##set B to all if missing
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
    #get j for heuristic and bound
    getR<-any(submat>=1:ncol(sCat))
    Hj<-ifelse(getR,1,0)  
    Bj<-Hj    }  
  
  #a subset of rows selected and some not changed
  if(type=="partial"){
    #subscripts where no action required
    Bfix<-c(B, rep(2,nrow(sCat)-length(B) ) )
    
    #create a new matrix by fixing/removing rows based on B
    submat<-rbind(sCat[which(B==1),], sCat[which(Bfix==2),])
    
    if(sum(Bfix==2)>1) submats<-rbind(sCat[which(B==1),], apply(sCat[which(Bfix==2),], 2, sort))
    if(sum(Bfix==2)==1) submats<-rbind(sCat[which(B==1),], sCat[which(Bfix==2),])
    
    #calculate cumsum over columns for customized mat
    cumcat.p<-apply(submat, 2, cumsum)
    
    #get j for heuristic
    Hj<-findj(cumcat.p)
    
    #calculate cumsum over columns for customized mat
    cumcat.s<-apply(apply(submat, 2, sort), 2, cumsum) 
    #get j for bound
    Bj<-findj(cumcat.s)
    
  }
  
  #all rows are indexed in B
  if(type=="full"){
    
    #create a new matrix by removing rows based on B
    submat<-sCat[which(B==1),]

    #calculate cumsum over columns for customized mat
    cumcat.p<-apply(submat, 2, cumsum)
    #get j for heuristic
    Hj<-findj(cumcat.p)
    
    #calculate cumsum over columns for customized mat
    cumcat.s<-apply(apply(submat, 2, sort), 2, cumsum) 
    #get j for bound
    Bj<-findj(cumcat.s)  
  }
  
  ##calculate the heuristic for all cases
  Hu=nrow(sCat)-(Hj-1)
  
  ##calculate the bound for all cases
  Bo=nrow(sCat)-(Bj-1)
  
  ###return the bound for true discoveries
  return(list("heuristic"=Hu, "bound"=Bo))
}
