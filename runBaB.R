#' @title Branch and bound algorithm implementation 
#'
#' @description performs B&B
#'
#' @param sCat  Category matrix, output of getCat function
#'  
#' @return A list of two objects, the heuristic and the bound
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
    step=1 
    
    ##current solution
    cFD<-ssh
    
    ##tracking parent bound
    PBlist<-c(ssb,ssb)
    
    ##initial branched
    branches<-list(1,0)
    
    outlist<-list()
    #search with B&B  until maximum step 
    while(step<=nMax & length(branches)>0){
      
      B<-branches[[1]]
      
      #run single step on the branch
      SSloc<-singleStep(sCat,B)
      outlist[[step]]<-c(SSloc$heuristic, SSloc$bound, B)
      
      #get the bound for branch
      Cbound<-min(PBlist[1], SSloc$bound)
      
      #remove the checked branch
      branches[[1]]<-NULL
      #remove the corresponding parent bound for checked branch
      PBlist<-PBlist[-1]
      
      #branching should continue 
      if( Cbound > cFD ){
        ##update the heuristic if a better one is found
        cFD<-max(cFD,SSloc$heuristic)
        
        #cont. branching if possible, add new branched to list
        if(length(B)<nrow(sCat) ){
          #avoid repeating main branch
          if( !(sum(B)==nrow(sCat)-1) ){
            branches[[length(branches)+1]]<-c(B,1)
            PBlist[length(PBlist)+1]<-Cbound
          }
          
          #avoid branch with all zero
          if( !(sum(B)==0 & length(B)==(nrow(sCat)-1)) ){
            branches[[length(branches)+1]]<-c(B,0) 
            PBlist[length(PBlist)+1]<-Cbound
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
    
    mBd<-max(cFD, max(BdinQ))
    
    #return solution
    return(list("FD"=cFD,"Bd"=mBd, "Step"=step, "nQ"=length(branches), "allout"=outlist)) 
}

