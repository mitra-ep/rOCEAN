#' @title OCEAN algorithm 
#'
#' @description Calculates heuristic and lower bound for the true discovery proportion (TDP) in 3 scales for
#' a specified two-way feature set (Algorithm 1 in the reference).
#' The input is either two omics data sub-matrices or the pre-calculated matrix of p-values for pairwise associations.
#' In case the result is not exact, the function adopts branch and bound (Algorithm 2 in the reference), if \code{nMax} allows.
#'
#' @param om1,om2  Matrix; Subsets of two omics data sets where rows are the features and columns are samples.
#' The rows of the two matrices would define the two-way feature set of interest.
#' 
#' @param gCT Vector; Parameters of the global closed testing, output of simesCT function.
#' 
#' @param scale Optional character vector; It specifies the scale of TDP quantification.
#' Possible choices are "pair" (pair-TDP), "col" (col-TDP ) and "row" (for row-TDP').
#' If not specified, all three scales are returned.
#' 
#' @param mps Optional matrix of p-values; A sub-matrix of pairwise associations, representing the two-way feature set of interest.
#' If provided, \code{om1} and \code{om2} are not required. If not provided, matrix of pairwise associations will be
#' derived from \code{om1} and \code{om2} based on Pearson's correlation.
#' 
#' @param nMax Maximum number of steps for branch and bound algorithm, if set to 1 branch and bound
#' is skipped even if the result is not exact. The default value is a 100. The algorithm may
#' stop before the \code{nMax} is reached if it converges sooner.
#' 
#' @param verbose Logical; if \code{TRUE}, progress messages will be displayed during the function's execution. Default is \code{TRUE}.
#' 
#' @return TDP is returned for the specified scales, along with number of steps taken and
#' convergence status for branch and bound algorithm.
#'
#' @author Mitra Ebrahimpoor
#'
#' \email{m.ebrahimpoor@@lumc.nl}
#' 
#' @seealso \link{simesCT}
#'  \link{pairTDP}
#'  \link{runbab}
#'
#' @export
#' 
#' 

ocean<-function(om1, om2, gCT, scale=c("pair", "row", "col"),
                  mps, nMax=100, verbose = TRUE) {

  #initiate
  concp=gCT[2]
  
  if(missing(scale)) {
    scale=c("pair", "row", "col")
  }
  #conditional print messages
  cat_if_verbose <- function(...) {
    if (verbose) cat(...)
  }
  
  ##############################pairTDP
  ptdp<-NA
  
  if(any("pair" %in% scale) & missing(mps)) {
    #calculate vector of p values
    pps<-corPs(om1, om2, type="Vec", pthresh=concp)
    ptdp<-pairTDP(pps, nrow(om1) * nrow(om2), gCT)
    cat_if_verbose("pair-TDP done. \n")
  }
  
  if(any("pair" %in% scale) & !missing(mps)) {
    pps<-as.vector(mps[])
    ptdp<-pairTDP(pps, nrow(mps) * ncol(mps), gCT)
    cat_if_verbose("pair-TDP done. \n")
  }
  
    #calculate cormat if necessary
  if(any(c("row", "col") %in% scale) & missing(mps)) {
    # calculate matrix of p values
    cat_if_verbose("Generating correlation matrix... \n")
    mps<-corPs(om1, om2, type="Mat")
    cat_if_verbose("Correlation matrix ready. \n")
  }
  
  ##############################row-tdp
  if("row" %in% scale) {
    sCatr<-getCat(mps, gCT, scale="row")
    gc()
    cat_if_verbose("p-categories matrix for rows ready. \n")
    #run single step algorithm
    ssr<-singleStep(sCatr)
    #td if conclusive result
    rtdp<-ssr$bound / nrow(sCatr)
    stepr<-1
    
    #td if inconclusive and run BaB
    gc()
    if (ssr$bound != ssr$heuristic & nMax > 1) {
      cat_if_verbose("Running BaB for row-TDP... \n")
      bboutr<-runbab(sCatr, ssr$heuristic, ssr$bound, nMax=nMax)
      stepr<-bboutr$Step
      bboutr$H<-bboutr$H / nrow(sCatr)
      bboutr$B<-bboutr$B / nrow(sCatr)
    }
    cat_if_verbose("row-TDP done. \n")
  } else {
    rtdp<-NA
    stepr<-NA
  }
  
  ##############################col-tdp
  if("col" %in% scale) {
    sCatc<-getCat(mps, gCT, scale="col")
    gc()
    cat_if_verbose("p-categories matrix for columns ready. \n")
    #run single step algorithm
    ssc<-singleStep(sCatc)
    
    #td if conclusive result
    ctdp<-ssc$bound / nrow(sCatc)
    stepc<-1
    
    #td if inconclusive and run BaB
    gc()
    if(ssc$bound != ssc$heuristic & nMax > 1) {
      cat_if_verbose("Running BaB for column-TDP... \n")
      bboutc<-runbab(sCatc, ssc$heuristic, ssc$bound, nMax=nMax)
      stepc<-bboutc$Step
      bboutc$H<-bboutc$H / nrow(sCatc)
      bboutc$B<-bboutc$B / nrow(sCatc)
    }
    cat_if_verbose("column-TDP done. \n")
  } else {
    ctdp<-NA
    stepc<-NA
  }
  
  #arrange items to return
  #result for row
  outr<-if(!is.na(stepr) & stepr>1) {
    c("rHeuristic"=bboutr$H, "rBound"=bboutr$B, "nStep"=stepr)
  } else if(!is.na(stepr) & stepr==1) {
    c("row-TDP"=rtdp, "nStep"=stepr)
  } else {
    NA
  }
  
  #result for col
  outc<-if(!is.na(stepc) & stepc>1) {
    c("cHeuristic"=bboutc$H, "cBound"=bboutc$B, "nStep"=stepc)
  } else if(!is.na(stepc) & stepc==1) {
    c("column-TDP"=ctdp, "nStep"=stepc)
  } else {
    NA
  }
  
  out<-list("Pairs"=c("pTDP"=ptdp),
              "Rows"=outr,
              "Columns"=outc)
  
  #remove NAs from item to return
  out<-out[!is.na(out)]
  
  return(out)
}
