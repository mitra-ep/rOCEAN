#' @title Calculates categories for a selected ROI 
#'
#' @description Calculates p-value categories from row data
#'
#' @param ps Matrix of p-values for the region of interest, representing association between 
#' two omics pathways
#' 
#' @param gCT Parameters of the global closed testing provided as the output of simesCT function
#'
#' @param m Size of the original p-value matrix (number of rows × number of columns) before 
#' selection of the region based on pathways.
#' 
#' @param scale Scale of the quantification, a character string. Possible choices are "col" and "row".
#' 
#' @return A matrix of categories where rows corresponds to the size of intersection
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
#' 

getCat<-function(ps, gCT, scale=c("col","row")){
  
  #parameters
  grandH=gCT[1]
  z=gCT[3]
  alpha=gCT[4]
  m=gCT[5]
  
  ##apply inversion if required
  if(scale=="row") sq.cat<-ceiling(ps[]*grandH/alpha)
  if(scale=="col") {
    ps<-t(ps)
    sq.cat<-ceiling(ps[]*grandH/alpha)}
  
  #get the size of categories from concentrations
  catSize<-min(as.numeric(nrow(sq.cat))*as.numeric(nrow(sq.cat)),
               max(unlist(sq.cat), na.rm=T), z-m+grandH+1, z)
  
  #calculate cumulative num in categories
  if(catSize>=2){
  sCat<-t(apply(sq.cat, 1, function(x) {
    out<-numeric(catSize)
    for (i in x)
      if (0 < i & i <= catSize & !is.na(i)) out[i]<-out[i]+1
    cumsum(out)
  }))
  
  #sort rows
  ordcat<-order(apply(sCat, 1, function(row) max(row / 1:length(row))))
  sCat<-sCat[ordcat,]
    }
  
  #all zero category mat for cases where no signal is detected
  if(catSize<2) sCat<-matrix(data=0,nrow=as.numeric(nrow(ps)), ncol=catSize+1)
  
  #message when done
  cat("Categories matrix created. \n")
  
  #remove large objects
  gc()
  
  return(sCat)}
