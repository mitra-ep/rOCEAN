#' @title ocean algorithm 
#'
#' @description Estimates TDP in 3 scales for a given region of interest
#'
#' @param om1,om2  The omics datasets where rows are probes and columns are samples.
#'
#' @param p1,p2  Optional argument to select pathways based on indexes (or row-names) of the omics datasets
#' which define the pathways of interest corresponding to each omic dataset
#'
#' @param gCT Parameters of the global closed testing provided as the output of simesCT function 
#' 
#' @param scale Optional, will specify the scale of the quantification, a character string. Possible choices are "pair", "col" and "row".
#' If not provided, all scales are returned.
#' 
#' @param saveP If true, saves the correlation matrix
#' @return For pair, fd is returned. For row and col, FD is returned if single step is conclusive 
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

ocean<-function(om1, om2, gCT, scale=c("pair","row","col"),
                mps, cor=FALSE, nMax=100){
  
  #check the pvalue matrix
  if(cor==TRUE & missing(mps)) stop("Matrix of p-values not provided.")

  #parameters
  concp=gCT[2]

  if(missing(scale)){
    scale=c("pair","row","col") }
  
  ptdp<-NA
  #TD over pairs
  if("pair" %in% scale & missing(mps)){
    #calculate vector of p values
    pps<-corPs(om1, om2, type="Vec", pthresh=concp )
    ptdp<-pairTD(pps, nrow(om1)*nrow(om2), gCT)
   }
    cat("check 1")
  #run pairwise algorithm
  if("pair" %in% scale & !missing(mps)){
    pps<-as.vector(mps[])
    ptdp<-pairTD(pps, nrow(mps)*ncol(mps), gCT)
    }


  
  if(sum(c("row","col") %in% scale)>=1 & missing(mps)){
    #calculate matrix of p values
    mps<-corPs(om1, om2, type="Mat")
  }
  cat("corP done \n")
  
  if("row" %in% scale){
      sCatr<-getCat(mps, gCT, scale="row")
      gc()
      cat("scat ready \n")
      #run single step algorithm
      ssr<-singleStep(sCatr)
      cat("Row done \n")
      #td if conclusive result                    
      rtdp=1-(ssr$bound/nrow(sCatr) )
      stepr=1

      #td if inconclusive and run BaB   
      gc()
      if(ssr$bound!=ssr$heuristic){ 
          bboutr<-runbab(sCatr,ssr$heuristic,ssr$bound,nMax=nMax)
          stepr=bboutr$Step
          bboutr$sL<-1-(bboutr$sL/nrow(sCatr) )
          bboutr$sU<-1-(bboutr$sU/nrow(sCatr) ) 
          cat("Row BB done \n")}
     
      
        } else{ rtdp=NA
                stepr=NA
                
         }
  
  
    if("col" %in% scale){
        sCatc<-getCat(mps, gCT,scale="col")
        gc()
        #run single step algorithm
        ssc<-singleStep(sCatc)
        
        #td if conclusive result                    
        ctdp=1-(ssc$bound/nrow(sCatc) )  
        stepc=1
        cat("Col done \n")
        #td if inconclusive and run BaB  
        gc()
        if(ssc$bound!=ssc$heuristic){ 
            bboutc<-runbab(sCatc,ssc$heuristic,ssc$bound,nMax=nMax)
            stepc=bboutc$Step 
            bboutc$sL<-1-(bboutc$sL/nrow(sCatc) )
            bboutc$sU<-1-(bboutc$sU/nrow(sCatc) )
            cat("Col BB done \n")}
        
          }  else{ ctdp=NA
                   stepc=NA
                   }

    #arrange items to return
    #result for row
    if(!is.na(stepr) & stepr>1){
      outr<-c("rHeuristic"=bboutr$sL,"rBound"=bboutr$sU, "nStep"=stepr)
    }
    if(!is.na(stepr) & stepr==1){
      outr<-c("rTDP"=rtdp, "nStep"=stepr)}
    if(is.na(stepr)) outr<-NA
  
    #result for col
    if(!is.na(stepc) & stepc>1){
      outc<-c("cHeuristic"=bboutc$sL,"cBound"=bboutc$sU, "nStep"=stepc)
    }
    if(!is.na(stepc) & stepc==1){
      outc<-c("cTDP"=ctdp, "nStep"=stepc)}
    if(is.na(stepc)) outc<-NA
  
    out<-list("Pairs"=c("pTDP"=ptdp),
              "Rows"=outr,
              "Columns"=outc  )
      
    ###remove NAs from item to return
    out<-out[!is.na(out)]

    return(out)
  }
