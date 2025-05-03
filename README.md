
# OCEAN

OCEAN is a flexible feature set testing method for analysis of multi-omics. For a pair of omics, either matrix of pairwise associations is derived from a pair of pre-processed omics data (using corPs, embedded in ocean) or the p-value matrix provided by the user is used.
This matrix is the input for simesCT, which calculates the closed testing parameters based on Simes local tests for a given alpha. Then for any given two-way feature set, 3 error rates are calculates using ocean function. 
* TDP: proportion of omics pairs that are associated.
* row-TDP: proportion of rows that are associated with at least one column.
* column-TPD: proportion of columns that are associated with at least one row.
TDP is calculated as an extension of SEA algorithm (see SEA R-package). For row-TDP and column-TDP a lower-bound (B) and heuristic (H) is calculated. Optionally, if the results are not exact (B not equal to H), it is possible to run branch and bound algorithm (using runBaB, embedded in ocean) to get an exact result. There are no limits on the number of feature sets being tested and the family-wise error rate is always controlled at level aplpha as set in the first step.

# Installation

You can install the development version of OCEAN from
[GitHub](https://github.com/) with:

``` r
install.packages("devtools")

#install the package from GitHub
devtools::install_github("mitra-ep/rOCEAN")
```
Also, a CRAN version is currently availble, which can be installed using:

``` r
install.packages("rOCEAN")

```
# Simulate data

To help you get started, here’s how to simulate a very simple and small dataset of p-values that you can use to test the functions of package.
This code generates a 1200 x 1000 matrix of p-values, with some signal intentionally added.

``` r
#number of feature per omic data set
n_cols<-1000
n_rows<-1200

#'#random matrix of p-values
set.seed(1258)
pvalmat<-matrix(runif(n_rows*n_cols, min=0, max=1)^3, nrow=n_rows, ncol=n_cols)
```

# CT parameters

Step one is to calculate the closed testing parameters:
```{r}
library(rOCEAN)
gCT<-simesCT(mps=pvalmat, m=nrow(pvalmat)*ncol(pvalmat))
```
The simesCT() function calculates five key parameters needed for downstream TDP estimation. For more details on these parameters, see:
Meijer, R. J., & Goeman, J. J. (2019). *Hommel’s procedure in linear time*. _Biometrika_, **106**(2), 483–489. [https://doi.org/10.1093/biomet/asz006](https://doi.org/10.1093/biomet/asz006)

These parameters are independent of any specific feature set and only need to be computed once for a given dataset and $\alpha$ level. You can reuse gCT for multiple runs of ocean() with different feature sets.

Note: Depending on the size of your data, this step may take time. For large datasets, consider applying a p-value threshold or using pre-computed p-value matrix (mps) to save time.


# TDP calculation

Now using the simulated p-value matrix, you can test the core functionality of rOCEAN based on a small two-way feature set.
Foolowing code estimates the true discovery proportion (TDP) at three levels, for the given feature set:

```{r}
library(rOCEAN)
#define a two-way feature set
subpmat <- pvalmat[1:40, 10:75]

#apply ocean function
out <- ocean(mps = subpmat, gCT = gCT, nMax = 2)
```
This is the output.

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

In the example above `nMax=2` so only 2 steps of BaB were applied for column-TDP, you can increase nMax for tighter bounds or leave it at the default of 100 or much higher for full refinement.

```{r}
#apply ocean function
out <- ocean(mps = subpmat, gCT = gCT, nMax = 100)
```
    #> p-categories matrix for columns ready. 
    #> Running BaB for column-TDP... 
    #> column-TDP done.
    #> $Columns
    #>  cHeuristic      cBound       nStep 
    #>   0.3839286   0.3750000 100.0000000


