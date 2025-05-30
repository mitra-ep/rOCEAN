#' @title Calculate cumulative p-categories
#'
#' @description Calculates cumulative p-categories for a given matrix of p-values.
#'
#' @param mps Matrix of p-values, representing pairwise associations between two feature sets.
#' 
#' @param gCT Parameters of the global closed testing, which is the output of simesCT function.
#'
#' @param scale Scale of the quantification, a character string. Possible choices are "col" and "row".
#' 
#' @return Matrix of p-categories.
#'
#' @seealso \link{simesCT}
#'
#' @examples
#' 
#' #number of features per omic data set
#' n_cols<-100
#' n_rows<-120
#' 
#' #random matrix of p-values
#' set.seed(1258)
#' pvalmat<-matrix(runif(n_rows*n_cols, min=0, max=1)^5, nrow=n_rows, ncol=n_cols)
#' 
#' #calculate CT parameters
#' gCT<-simesCT(mps=pvalmat, m=nrow(pvalmat)*ncol(pvalmat))
#' 
#' #define the two-way feature set
#' subpmat<-pvalmat[61:75,81:100]
#' 
#' #calculate p-categories matrix for feature set by rows
#' rCat<-getCat(mps=subpmat, gCT, scale="row")
#'
#' #calculate p-categories matrix for feature set by columns
#' cCat<-getCat(mps=subpmat, gCT, scale="col")
#'
#' @export
#' 
#' 

getCat<-function(mps, gCT, scale=c("col","row")){
  
  #donotrun if no signal
  if(gCT[3]==0) stop("No discoveries in this data!")
  
  #parameters
  grandH=gCT[1] 
  z=gCT[3]
  alpha=gCT[4]
  m=gCT[5]
  
  
  scale<-match.arg(scale)
  
  ##apply inversion if required
  if(scale=="row") sq.cat<-ceiling(mps[]*grandH/alpha)
  if(scale=="col") {
    mps<-t(mps)
    sq.cat<-ceiling(mps[]*grandH/alpha)
    }
  
  #get the size of categories from concentrations
  catSize<-min(as.numeric(nrow(sq.cat))*as.numeric(nrow(sq.cat)),
               max(unlist(sq.cat), na.rm=T), z-m+grandH+1, z)

  #calculate cumulative num in categories
  if(catSize>=2){
    sCat <- t(apply(sq.cat, 1, function(x) {
      out <- numeric(catSize)
      for (i in x) {
        if ((0 < i) && (i <= catSize) && !is.na(i)){
          out[i] <- out[i] + 1
        }
      }
      cumsum(out)
    }))
    
  
  #sort rows
  ordcat<-order(apply(sCat, 1, function(row) max(row / 1:length(row))))
  sCat<-sCat[ordcat,]
    }

  #all zero category mat for cases where no signal is detected
  if(catSize<2) sCat<-matrix(data=0,nrow=as.numeric(nrow(mps)), ncol=catSize+1)

  #remove large objects
  gc()
  
  return(sCat)}
