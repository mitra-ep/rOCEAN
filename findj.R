#' @title Find J 
#'
#' @description For a matrix with elements presented as v_jk,
#' searches for the smallest j where v{jk>=k 
#'
#' @param mat input matrix
#'  
#' @return value of J
#'
#' @author Mitra Ebrahimpoor
#'
#' \email{m.ebrahimpoor@@lumc.nl}
#'
#' 

findj <- function(mat){
  r<- nrow(mat)
  l<- ncol(mat)
  
  #if no mat_ij>=k
  res=r
  
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