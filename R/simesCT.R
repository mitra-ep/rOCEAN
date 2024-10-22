#' @title Calculates parameters of global closed testing based on Simes test 
#'
#' @description Calculates five parameters from closed testing with Simes local tests based on row data.
#'  Grand H: size of the largest set which is not rejected by closed testing and z: size of concentration set.
#'  These values are unique per dataset/alpha-level combination and do not depend on pathways. Calculation may
#'  be somewhat long depending on the size of datasets and PC configurations.
#'
#' @param om1,om2 Two omics datasets where rows are probs and columns are samples.
#' 
#' @param mps,m Optional, pre-calculated matrix/vector of pairwise associations and the size.
#' To save time in calculation of parameters, one can apply a threshold such as the type I error to rewove larger p-values.
#' If a threshold is used size of matrix and m will not match. m should always be the size of the original matrix of associations.
#' 
#' @param alpha type I error rate
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
  concP=sp[z]  
  #remove large objects
  gc()

  return(c("grandH"=grandH,"concP"=concP,"z"=z,"alpha"=alpha, "m"=m) )}
