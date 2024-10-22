#' @title Calculates p-value  
#'
#' @description Calculates pairwise matrix of p-values based on Pearson's correlation test for two matrices.
#' To gain speed, the matrix is split into several smaller blocks.
#' 
#' @param om1,om2  Subsets of omics data where rows are the probes and columns are samples.
#' The rows of the two matrices should define the two-way feature set of interest.
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
#' @importFrom ff ff
#' 
#' @export
#' 

corPs<-function(om1, om2,type=c("Mat","Vec"), pthresh){
   
  #check
  if(ncol(om1)!=ncol(om2))
    stop("The two matrices must have the same number of columns (i.e. same number of subjects).")
 
  #pval from pearson corr
  cor2p<-function(r,n=ncol(om1)){
    t<-(r*sqrt(n-2))/sqrt(1-r^2)
    p<-2*(1 - pt(abs(t),(n-2)))
    return(p)
  }
  
  #pearson corr while avoiding zero-var
  custom_cor_mat<-function(mat1, mat2) {
    #columns with zero-var
    var1<-apply(mat1, 2, var)
    var2<-apply(mat2, 2, var)
    
    non_zero_var1<-which(var1 != 0)
    non_zero_var2<-which(var2 != 0)
    
    mat1_filtered<-mat1[, non_zero_var1, drop = FALSE]
    mat2_filtered<-mat2[, non_zero_var2, drop = FALSE]
    
    #scale matrices
    mat1_scaled<-scale(mat1_filtered)
    mat2_scaled<-scale(mat2_filtered)
    
    #get pearson corr
    cor_matrix<-cor(mat1_scaled, mat2_scaled)
    
    #full result matrix with 0
    result<-matrix(0, ncol(mat1), ncol(mat2))
    
    # fill in the final mat
    result[non_zero_var1, non_zero_var2]<-cor_matrix
    
    #allocate names
    rownames(result)<-colnames(mat1)
    colnames(result)<-colnames(mat2)
    return(result)
  }

    #find factors function
    mfact <- function(n){
      one.to.n <- seq_len(n)
      fact<-one.to.n[(n %% one.to.n) == 0]
      max(fact[fact<30])
    }
    

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
      COR<-custom_cor_mat(t(om1[G1,]), t(om2[G2,]) )
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
        COR<-custom_cor_mat(t(om1[G1,]), t(om2[G2,]) )
        CORP<-cor2p(COR)
        corOUT[G1, G2]<-CORP
        COR<-NULL
        gc()
      }
      
      
    }

  return(corOUT)
}
