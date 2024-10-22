#' @title OCEAN algorithm 
#'
#' @description Estimates TDP bounds in 3 scales for a given pathway.
#' In case the results is not exact, the function adopts branch and bound if nMax allows.
#'
#' @param om1,om2  Subsets of omics data where rows are the probes and columns are samples.
#' The rows of the two matrices should define the two-way feature set of interest.
#' 
#' @param gCT Parameters of the global closed testing provided as the output of simesCT function 
#' 
#' @param scale Optional, will specify the scale of the quantification, a character string. Possible choices are "pair", "col" and "row".
#' If not provided, all scales are returned.
#' 
#' @param mps Optional, a sub-matrix of the matrix of pairwise associations. If provided, om1 and om2 are not required.
#' If not provided, matrix of pairwise associations will be derived from om1 and om2 based on Pearson's correlation.
#' 
#' @param nMax Maximum number of steps for branch and bound, if set to 1 branch and bound
#' is skipped even if the result is not exact. 
#' 
#' @return TDP is returned for the scales defined by scale. 
#'
#' @author Mitra Ebrahimpoor
#'
#' \email{m.ebrahimpoor@@lumc.nl}
#'
#' @export
#' 
#' 

ocean <- function(om1, om2, gCT, scale = c("pair", "row", "col"),
                  mps, nMax = 100) {
  
  # check the pvalue matrix
  if (missing(om1) & missing(om2)){
    if (missing(mps)) stop("Matrix of p-values not provided.")
     }
    
  
  # parameters
  concp = gCT[2]
  
  if (missing(scale)) {
    scale = c("pair", "row", "col")
  }
  
  ptdp <- NA
  # TD over pairs
  if (any("pair" %in% scale) & missing(mps)) {
    # calculate vector of p values
    pps <- corPs(om1, om2, type = "Vec", pthresh = concp)
    ptdp <- pairTD(pps, nrow(om1) * nrow(om2), gCT)
  }
  
  # run pairwise algorithm
  if (any("pair" %in% scale) & !missing(mps)) {
    pps <- as.vector(mps[])
    ptdp <- pairTD(pps, nrow(mps) * ncol(mps), gCT)
  }
  
  if (any(c("row", "col") %in% scale) & missing(mps)) {
    # calculate matrix of p values
    mps <- corPs(om1, om2, type = "Mat")
  }
  
  cat("corP done \n")
  
  if ("row" %in% scale) {
    sCatr <- getCat(mps, gCT, scale = "row")
    gc()
    cat("scat ready \n")
    # run single step algorithm
    ssr <- singleStep(sCatr)
    cat("Row done \n")
    # td if conclusive result
    rtdp <- ssr$bound / nrow(sCatr)
    stepr <- 1
    
    # td if inconclusive and run BaB
    gc()
    if (ssr$bound != ssr$heuristic & nMax > 1) {
      bboutr <- runbab(sCatr, ssr$heuristic, ssr$bound, nMax = nMax)
      stepr <- bboutr$Step
      bboutr$H <- bboutr$H / nrow(sCatr)
      bboutr$B <- bboutr$B / nrow(sCatr)
      cat("Row BB done \n")
    }
    
  } else {
    rtdp <- NA
    stepr <- NA
  }
  
  if ("col" %in% scale) {
    sCatc <- getCat(mps, gCT, scale = "col")
    gc()
    # run single step algorithm
    ssc <- singleStep(sCatc)
    
    # td if conclusive result
    ctdp <- ssc$bound / nrow(sCatc)
    stepc <- 1
    cat("Col done \n")
    # td if inconclusive and run BaB
    gc()
    if (ssc$bound != ssc$heuristic & nMax > 1) {
      bboutc <- runbab(sCatc, ssc$heuristic, ssc$bound, nMax = nMax)
      stepc <- bboutc$Step
      bboutc$H <- bboutc$H / nrow(sCatc)
      bboutc$B <- bboutc$B / nrow(sCatc)
      cat("Col BB done \n")
    }
    
  } else {
    ctdp <- NA
    stepc <- NA
  }
  
  # arrange items to return
  # result for row
  if (!is.na(stepr) & stepr > 1) {
    outr <- c("rHeuristic" = bboutr$H, "rBound" = bboutr$B, "nStep" = stepr)
  }
  if (!is.na(stepr) & stepr == 1) {
    outr <- c("rTDP" = rtdp, "nStep" = stepr)
  }
  if (is.na(stepr)) outr <- NA
  
  # result for col
  if (!is.na(stepc) & stepc > 1) {
    outc <- c("cHeuristic" = bboutc$H, "cBound" = bboutc$B, "nStep" = stepc)
  }
  if (!is.na(stepc) & stepc == 1) {
    outc <- c("cTDP" = ctdp, "nStep" = stepc)
  }
  if (is.na(stepc)) outc <- NA
  
  out <- list("Pairs" = c("pTDP" = ptdp),
              "Rows" = outr,
              "Columns" = outc)
  
  # remove NAs from item to return
  out <- out[!is.na(out)]
  
  return(out)
}
