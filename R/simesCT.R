#' @title Closed testing with Simes
#'
#' @description Calculates five parameters from closed testing with Simes local tests based on raw data.
#'  These parameter are unique per data/alpha-level combination and do not depend on feature sets. Calculation may
#'  be somewhat long depending on the size of data sets and PC configurations.
#'
#' @param om1,om2 Two omics data sets where rows are features and columns are samples.
#' 
#' @param mps,m Optional, pre-calculated matrix/vector of pairwise associations and the size.
#' To save time in calculation of parameters, a threshold such as the type I error may be applies to remove larger p-values.
#' If a threshold is used, size of matrix and \code{m} will not match. \code{m} should always be the size of the matrix of associations
#' (number of features in \code{om1}  X number of features in \code{om2}).
#' 
#' @param alpha type I error rate, default value is 0.05.
#'
#' @return Vector of integers: grand H value, concentration p-value, size of concentration set z,
#' size of the original pair-wise associations matrix and the type I error level used in calculations.
#'
#' @author Mitra Ebrahimpoor
#'
#' \email{m.ebrahimpoor@@lumc.nl}
#'
#' @references See more details in "Hommel's procedure in linear time" doi: 10.1002/bimj.201700316.
#'
#' @export
#' 

simesCT<-function(om1, om2, mps, m, alpha=0.05){
  
  ##checks
  if(!missing(om1) && !missing(om2)){
    if(ncol(om1)!=ncol(om2)) stop("ncol of the matrices should match!")
  }else{
      if (missing(mps)) stop("No data pair or matrix of p-values provided.")}

  #mps not provided
  if (missing(mps)){
    #get all pvals after filtering by alpha
  m<-as.numeric(nrow(om1))*as.numeric(nrow(om2))
  sp<-corPs(om1, om2, type= "Vec", pthresh=alpha)}
  else{
    sp<-as.vector(mps<alpha)
  }
  sp<-sort(sp[])
  k<-length(sp)
  grandH<- m-max(0, ceiling(max(1:k - (m-1:k) * sp / (alpha - sp))))

  #get size of concentration set
  z<-ifelse(grandH==m, 0, min(which(sp*grandH <= (1:k - m + grandH + 1) * alpha)))
  
  # if no signal in data, set concP to alpha as only used for filtering
  concP <- ifelse(z == 0, 0.05, sp[z])
  
  #remove large objects
  gc()

  return(c("grandH"=grandH,"concP"=concP,"z"=z,"alpha"=alpha, "m"=m) )}
