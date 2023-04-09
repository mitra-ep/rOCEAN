#' @title Branch and bound algorithm implementation 
#'
#' @description performs B&B
#'
#' @param sCat  Category matrix, output of getCat function
#'  
#' @return A list of two objects, the heuristic, the bound, number of steps
#' taken and convergence status
#'
#' @author Mitra Ebrahimpoor
#'
#' \email{m.ebrahimpoor@@lumc.nl}
#'
#' @seealso
#'
#' @references
#'
#' @examples
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
    while(step<nMax & length(branches)>0){
      
      B<-branches[[length(branches)]]
      
      #run single step on the branch
      SSloc<-singleStep(sCat,B)
    
      #get the bound for branch
      Cbound<-min(SSloc$bound,PBlist[length(branches)])
      
      #remove the checked branch
      branches[[length(branches)]]<-NULL
      
      #remove the corresponding parent bound for checked branch
      PBlist<-PBlist[-length(branches)]
      
      #branching should continue 
      if( Cbound > sL ){
        ##update the heuristic if a better one is found
        sL<-max(sL,SSloc$heuristic)
        
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
      
    }else BdinQ<-0
    
    sU<-max(sL, max(BdinQ))
    
    
    #return solution
    return(list("sl"=sL,"sU"=sU, "Step"=step, "Conv"=length(branches)==0) ) 
  }
  
