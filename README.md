
# OCEAN

OCEAN is a flexible feature set testing method for analysis of
multi-omics. For a pair of omics, either matrix of pairwise associations
is derived from a pair of pre-processed omics data (using corPs,
embedded in ocean) or the p-value matrix provided by the user is used.
This matrix is the input for simesCT, which calculates the closed
testing parameters based on Simes local tests for a given alpha. Then
for any given two-way feature set, 3 error rates are calculates using
ocean function. The 3 error rates are TDP: proportion of omics pairs
that are truly associated; row-TDP: proportion of rows that are
associated with at least one column; column-TPD: proportion of columns
that are associated with at least one row. TDP is calculated as an
extension of SEA algorithm (see SEA R-package). For row-TDP and
column-TDP a lower-bound (B) and heuristic (H) is calculated.
Optionally, if the results are not exact (B not equal to H), it is
possible to run branch and bound algorithm (using runBaB, embedded in
ocean) to get an exact result. There are no limits on the number of
feature sets being tested and the family-wise error rate is always
controlled at level aplpha as set in the first step.

# Installation

You can install the development version of OCEAN from
[GitHub](https://github.com/) with:

``` r
install.packages("devtools")

#install the package from GitHub
devtools::install_github("mitra-ep/rOCEAN")
```

# Simulate data

``` r
library(rOCEAN)

#number of feature per omic data set
n_cols<-1000
n_rows<-1200

#'#random matrix of p-values
set.seed(1258)
pvalmat<-matrix(runif(n_rows*n_cols, min=0, max=1)^3, nrow=n_rows, ncol=n_cols)
```

# CT parameters

Calculate the closed testing parameters:

# TDP calculation

Calculate TDPs for an imaginary two-way feature set.

    #> pair-TDP done. 
    #> p-categories matrix for rows ready. 
    #> row-TDP done. 
    #> p-categories matrix for columns ready. 
    #> Running BaB for column-TDP... 
    #> column-TDP done.
    #> $Pairs
    #>        pTDP 
    #> 0.007936508 
    #> 
    #> $Rows
    #>   row-TDP     nStep 
    #> 0.5061728 1.0000000 
    #> 
    #> $Columns
    #> cHeuristic     cBound      nStep 
    #>  0.3839286  0.3750000  2.0000000

In the example above `nMax=2` so only 2 steps of BaB were applied for
column-TDP, we can increase the number of steps to make the bound
narrower:

    #> p-categories matrix for columns ready. 
    #> Running BaB for column-TDP... 
    #> column-TDP done.
    #> $Columns
    #>  cHeuristic      cBound       nStep 
    #>   0.3839286   0.3750000 100.0000000

Calculate TDPs for a case where the initial outcome is unsure and branch
and bound will be adopted for row-TDP.

    #> p-categories matrix for rows ready. 
    #> Running BaB for row-TDP... 
    #> row-TDP done.
    #> $Rows
    #>  rHeuristic      rBound       nStep 
    #>   0.9600000   0.9466667 100.0000000
