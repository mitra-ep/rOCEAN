#' @title Branch and bound algorithm implementation 
#'
#' @description Performs B&B when the bound are not exact
#' 
#' @param sCat  Category matrix, output of getCat function
#'  
#' @param ssh current Heuristic as provided by SingleStep function
#'
#' @param ssb current Bound as provided by SingleStep function
#' 
#' @param nMax Maximum number of steps for the algorithm, the algorithm may stop sooner if it converges.
#'  
#' @return A list, including the heuristic and the bound for the number of true discoveries, along with number of steps
#' taken and convergence status.
#'
#' @seealso \link{getCat}
#'  \link{singleStep}
#' 
#' @examples
#' 
#' #number of features per omic data set
#' n_cols<-100
#' n_rows<-120
#' 
#' #random matrix of p-values
#' set.seed(1258)
#' pvalmat<-matrix(runif(n_rows*n_cols, min=0, max=1)^4, nrow=n_rows, ncol=n_cols)
#' 
#' #calculate CT parameters
#' gCT<-simesCT(mps=pvalmat, m=nrow(pvalmat)*ncol(pvalmat))
#' 
#' #define the two-way feature set
#' subpmat<-pvalmat[1:10,31:40]
#' 
#' #calculate p-categories matrix for feature set by rows
#' rCat<-getCat(mps=subpmat, gCT, scale="row")
#' 
#' #calculate the heuristic and bound
#' SSout<-singleStep(rCat)
#' 
#' #run branch nd bound
#' runbab(rCat, SSout$heuristic, SSout$bound, nMax=800)
#' 
#' @export
#' 
#' 

runbab<-function(sCat, ssh, ssb, nMax=100){

    #initial arguments
    ##step count
    step=0 
    
    ##current solution
    sL<-ssh
    
    ##tracking parent bound
    PBlist<-c(ssb,ssb)
    
    ##initial branched
    branches<-list(1,0)
  
    #search with B&B  until maximum step 
    while(step<nMax){
     #stop the loop when converged
       if(length(branches) == 0) {
        break
      }
      B<-branches[[length(branches)]]
      
      #run single step on the branch
      SSloc<-singleStep(sCat,B)
    
      #get the bound for branch
      Cbound<-max(SSloc$bound,PBlist[length(branches)])
      
      #remove the checked branch
      branches[[length(branches)]]<-NULL
      
      #remove the corresponding parent bound for checked branch
      PBlist<-PBlist[-length(PBlist)]
      
      #branching should continue 
      if( Cbound < sL ){
        ##update the heuristic if a better one is found
        sL<-min(sL,SSloc$heuristic)
        
        #cont. branching if possible, add new branched to list
        if(length(B)<nrow(sCat) ){
          #avoid repeating main branch
          if( !(sum(B)==nrow(sCat)-1) ){
            branches[[length(branches)+1]]<-c(B,1)
            PBlist[length(PBlist)+1]<-SSloc$bound
          }
          
          #avoid branch with all zero
          if( !(sum(B)==0 & length(B)==(nrow(sCat)-1)) ){
            branches[[length(branches)+1]]<-c(B,0) 
            PBlist[length(PBlist)+1]<-SSloc$bound
          }
        }
        
      }
      #update step count
      step<-step+1
    }
    
    #update the estimate for bound based on branches left in queue
    if(length(branches)>0){
      #single step
      BdinQ<-unlist(lapply(branches, function(b) singleStep(sCat,b)$bound))
      #compare to parent bound
      BdinQ<-pmin(BdinQ,PBlist)
      
    }else BdinQ<-nrow(sCat)
    
    sU<-min(sL, min(BdinQ))
    
    
    #return solution
    return(list("B"=sU, "H"=sL, "Step"=step, "Conv"=length(branches)==0) ) 
  }
  
