#' @title Calculates p-value by spiliting 
#'
#' @description Calculate p-value matrix of Pearson correlation test for two matrices by splitting 
#' the rows to (25 or less) smaller blocks. If the size is larger than 1000X1000,
#' 
#' @param om1,om2  Two omics datasets where rows are probs and columns are samples.
#'
#' @param p1,p2  Optional argument to select pathways based on indexes (or row-names) of the omics datasets
#' which define the pathways of interest corresponding to each omic dataset
#' 
#' @param type  Two options are available. Mat: Calculate the correlation of subsets and return a
#' matrix; Vec: calculate the full correlation matrix, subset by the given threshold and return a
#' vector of p-values.
#' 
#' @param pthresh The threshold by which the p-values are filtered (p>pthresh is removed).
#' Default is 0.05
#' 
#' @return Either an ff_matrix or a vector, as determined by type parameter.
#'
#' @author Mitra Ebrahimpoor
#'
#' \email{m.ebrahimpoor@@lumc.nl}
#'
#' @seealso
#'
#'
#' @importFrom ff ff
#' 
#' 

corPs<-function(om1, om2, p1, p2,
                type=c("Mat","Vec"), pthresh){
   
    
    #check pathways
    if(!is.numeric(p1) & sum(p1 %in% rownames(om1))==0) stop("Cannot select from om1 with p1.")
    if(!is.numeric(p2) & sum(p2 %in% rownames(om2))==0) stop("Cannot select from om2 with p2.")
    
    if(is.numeric(p1) & any(p1>nrow(mo1)) ) stop("Cannot select from om1 with p1.")
    if(is.numeric(p2) & any(p2>nrow(mo2)) ) stop("Cannot select from om2 with p2.")
  
  #pval from pearson cor
    cor2p<-function(r,n=ncol(om1)){
      t<-(r*sqrt(n-2))/sqrt(1-r^2)
      p<-2*(1 - pt(abs(t),(n-2)))
      return(p)
    }
    
    #find factors function
    mfact <- function(n){
      one.to.n <- seq_len(n)
      fact<-one.to.n[(n %% one.to.n) == 0]
      max(fact[fact<30])
    }
    
    #make selection
    om1<-om1[p1,]
    om2<-om2[p2,]
    
    #choose number of splits
    s1<-mfact(nrow(om1))
    s2<-mfact(nrow(om2))
    
    #split column numbers into several block
    SPLIT1<-split(1:nrow(om1), rep(1:s1, each = nrow(om1)/s1))
    SPLIT2<-split(1:nrow(om2), rep(1:s2, each = nrow(om2)/s2)) 
    
    ## create all unique combinations of block
    COMBS<-expand.grid(1:length(SPLIT1), 1:length(SPLIT2))
    COMBS<-unique(COMBS)
    
    
    if(type=="Vec"){
  
    #make an empty vector
    corOUT<-ff::ff(vmode = "single", length=0)
    
    ## iterate through each combination
    ## calculate correlation 
    ## convert correlation to p-value
    ## store them in the vector 
    for (i in 1:nrow(COMBS)) {
      G1<-SPLIT1[[COMBS[i,1]]]
      G2<-SPLIT2[[COMBS[i,2]]]
      flush.console()
      COR<-cor(t(om1[G1,]), t(om2[G2,]) )
      CORP<-cor2p(COR)
      CORP<-CORP[CORP<pthresh]
      corOUT<-c(corOUT[],CORP)
      COR<-NULL
      gc()
    } 
    
    }
    
    if(type=="Mat"){
      
      ##make an empty matrix 
      corOUT<-ff::ff(vmode = "single", dim = c(nrow(om1), nrow(om2)))

      ## iterate through each combination
      ## calculate correlation 
      ## convert correlation to p-value
      ## store them in the vector      
      cat("Generating correlation matrix. \n")
      for (i in 1:nrow(COMBS)) {
        G1<-SPLIT1[[COMBS[i,1]]]
        G2<-SPLIT2[[COMBS[i,2]]]
        flush.console()
        COR<-cor(t(om1[G1,]), t(om2[G2,]) )
        CORP<-cor2p(COR)
        corOUT[G1, G2]<-CORP
        COR<-NULL
        gc()
      }
      
      
    }

  return(corOUT)
}
