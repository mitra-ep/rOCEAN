
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
devtools::install_github("mitra-ep/OCEA")
```

# Simulate data

This is a basic example which shows you how to solve a common problem.
For the sake of this example, we will simulate a correlation matrix.

``` r
library(OCEAN)

#reproduce
set.seed(123)

#dimensions
n_om1<-400
n_om2<-500
n_samp<-50

#two random data sets
om1<-matrix(rnorm(n_om1*n_samp, mean=0, sd=1), nrow=n_om1, ncol=n_samp)
om2<-matrix(rnorm(n_om2*n_samp, mean=0, sd=1), nrow=n_om2, ncol=n_samp)

#add signal
signal<-matrix(rnorm(50*n_samp, mean=82, sd=1), nrow=50, ncol=n_samp)
om1[1:50, ]<-om1[1:50, ]+signal
om2[21:70, ]<-om2[21:70, ]+signal
```

# CT parameters

Calculate the closed testing parameters:

# CT parameters

Calculate TDPs for an imaginary two-way feature set.

    #> pair-TDP done. 
    #> Generating correlation matrix... 
    #> Correlation matrix ready. 
    #> p-categories matrix for rows ready. 
    #> row-TDP done. 
    #> p-categories matrix for columns ready. 
    #> column-TDP done.
    #> $Pairs
    #> pTDP 
    #>    0 
    #> 
    #> $Rows
    #> row-TDP   nStep 
    #>       0       1 
    #> 
    #> $Columns
    #> column-TDP      nStep 
    #>          0          1
