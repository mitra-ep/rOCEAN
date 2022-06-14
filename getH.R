#' @title Calculates grand H
#'
#' @description Calculates grand H value (size f largt set not rejected by closed testing) from row data  
#'
#' @param om1,om2 omics datasets
#'
#' @return Integer, the grand H value.
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

getH<-function(om1, om2, alpha=0.05){
  
  ##checks
  if(ncol(om1)!=ncol(om2)) stop("ncol of the matrices should math!")

  ##calculate the pval mat for selection 
  sq<-getPs(om1, om2, p1=path1, p2=path2, type= "Mat")
  gc()
  
  ##calculation of grandH
  #get all pvals and filter by alpha
  m<-nrow(om1)*nrow(om2)
  sp<-getPs(om1, om2, type= "Vec", pthresh=alpha)
  sp<-sort(sp)
  k<-length(sp)
  grandH<- m-max(0, ceiling(max(1:k - (m-1:k) * sp / (alpha - sp))))
  
  #remove large objects
  gc()
  
  return(grandH)}
