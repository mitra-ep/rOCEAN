#' @title Calculates p-value categories from row data  
#'
#' @description Calculates p-value categories from row data
#'
#' @param pathlist  \code{subject}.
#'
#' @return A datafram including the name of pathways and corresponding adjusted p-values.
#'
#' @author Mitra Ebrahimpoor
#'
#' \email{m.ebrahimpoor@@lumc.nl}
#'
#' @seealso
#'
#' \code{\link{}}, \code{\link{}}, \code{\link{}}
#'
#' @references
#'
#' @examples
#'
#'
#' @export
#'
#' @importFrom ff ff
#' 
#' 

getCat<-function(om1, om2, path1, path2, alpha=0.05){
  
  ##checks
  if(ncol(om1)!=ncol(om2)) stop("ncol of the matrices should math!")
  if(missing(path1)|missing(path2)) stop("Pathway indexes are not defined.")
  
  ##calculate the pval mat for selection 
  sq<-getPs(om1, om2, p1=path1, p2=path2, type= "Mat")
  gc()
  
  ##calculation of grandH
  #get all pvals and filter by alpha
  m<-nrow(om1)*nrow(om2)
  sp<-getPs(om1, om2, type= "Vec", pthresh=alpha)
  k<-length(sp)
  grandH<- m-max(0, ceiling(max(1:k - (m-1:k) * sp / (alpha - sp))))
  
  ##calculate p-value categories
  sq.cat<-ceiling(sq[]*grandH/alpha)
  
  #get the size of categories from concentrations
  z<-ifelse(grandH==m, 0, min(which(sp*grandH <= (1:k - m + grandH + 1) * alpha)))
  catSize<-min(nrow(sq)*ncol(sq), max(unlist(sq.cat)), z-m+grandH+1)
  
  cat("Creating category marix. \n")
  #calculate cumulative num in categories
  sCat<-t(apply(sq.cat, 1, function(x) {
    out<-numeric(catSize)
    for (i in x)
      if (0 < i & i <= catSize) out[i]<-out[i]+1
    cumsum(out)
  }))
  
  #sort rows
  ordcat<-order(apply(sCat, 1, function(row) max(row / 1:length(row))))
  sCat<-sCat[ordcat,]

  #remove large objects
  
  return(sCat)}
