#' @title getCT calculates parameters of closed testing procedure 
#'
#' @description Calculates two parameters for closed testing with Simes local tests from row data.
#'  Grand H: size of the largest set which is not rejected by closed testing and z: size of concentration set.
#'  These values are unique per dataset/alpha-level combination and do not depend on pathways. Calculation may
#'  be somewhat long depending on the size of datasets and PC configurations.
#'
#' @param om1,om2 Two omics datasets where rows are probs and columns are samples.
#' 
#' @param alpha The type I error allowed
#'
#' @return Vector of 3 integers. The grand H and z value along with the corresponding type I error level.
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
#' 

simesCT<-function(om1, om2, alpha=0.05){
  
  ##checks
  if(ncol(om1)!=ncol(om2)) stop("ncol of the matrices should math!")

  ##calculation of grandH
  
  #get all pvals after filtering by alpha
  m<-as.numeric(nrow(om1))*as.numeric(nrow(om2))
  sp<-corPs(om1, om2, type= "Vec", pthresh=alpha)
  sp<-sort(sp[])
  k<-length(sp)
  grandH<- m-max(0, ceiling(max(1:k - (m-1:k) * sp / (alpha - sp))))

  #get size of concentration set
  z<-ifelse(grandH==m, 0, min(which(sp*grandH <= (1:k - m + grandH + 1) * alpha)))
  concP=sp[z]  
  #remove large objects
  gc()
  m=as.numeric(nrow(om1))*as.numeric(nrow(om2))
  
  return(c("grandH"=grandH,"concP"=concP,"z"=z,"alpha"=alpha, "m"=m) )}
