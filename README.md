
# OCEAN for multi-omics

OCEAN is a flexible feature set testing method for analysis of multi-omics. For a pair of omics, either matrix of pairwise associations is derived from a pair of pre-processed omics data (using corPs, embedded in ocean) or the p-value matrix provided by the user is used.
This matrix is the input for simesCT, which calculates the closed testing parameters based on Simes local tests for a given alpha. Then for any given two-way feature set, 3 error rates are calculates using ocean function. 
* TDP: proportion of omics pairs that are associated.
* row-TDP: proportion of rows that are associated with at least one column.
* column-TPD: proportion of columns that are associated with at least one row.
  
TDP is calculated as an extension of SEA algorithm (see rSEA R-package). For row-TDP and column-TDP a lower-bound (B) and heuristic (H) is calculated. Optionally, if the results are not exact (B not equal to H), it is possible to run branch and bound algorithm (using runBaB, embedded in ocean) to get an exact result. There are no limits on the number of feature sets being tested and the family-wise error rate is always controlled at level aplpha as set in the first step.

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
# Simulated data

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

## CT parameters

Step one is to calculate the closed testing parameters:

```{r}
library(rOCEAN)
gCT<-simesCT(mps=pvalmat, m=nrow(pvalmat)*ncol(pvalmat))
```
The simesCT() function calculates five key parameters needed for downstream TDP estimation. For more details on these parameters, see:
Meijer, R. J., & Goeman, J. J. (2019). *Hommel’s procedure in linear time*. _Biometrika_, **106**(2), 483–489. [https://doi.org/10.1093/biomet/asz006](https://doi.org/10.1093/biomet/asz006)

These parameters are independent of any specific feature set and only need to be computed once for a given dataset and $\alpha$ level. You can reuse gCT for multiple runs of ocean() with different feature sets.

## TDP calculation

Now using the simulated p-value matrix, you can test the core functionality of rOCEAN based on a small two-way feature set.
Foolowing code estimates the true discovery proportion (TDP) at three levels, for the given feature set:

```{r}
library(rOCEAN)
#define a two-way feature set
subpmat <- pvalmat[1:40, 10:75]

#apply ocean function
out <- ocean(mps = subpmat, gCT = gCT, nMax = 2)
```
This is an example output. 
<pre><code class="language-markdown"> ```
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
 ``` </code></pre>
 
In the example above `nMax=2` so only 2 steps of BaB were applied for column-TDP, you can increase nMax for tighter bounds or leave it at the default of 100 or much higher for full refinement.

```{r}
#apply ocean function
out <- ocean(mps = subpmat, gCT = gCT, nMax = 100)
```
The results:
<pre><code class="language-markdown"> ```
    #> p-categories matrix for columns ready. 
    #> Running BaB for column-TDP... 
    #> column-TDP done.
    #> $Columns
    #>  cHeuristic      cBound       nStep 
    #>   0.3839286   0.3750000 100.0000000
 ``` </code></pre>
 
# Raw data

In practice, you may sart with raw data. 
n the context of rOCEAN, the two omics datasets (omics1 and omics2) should have the following structure:

* Rows represent the features (e.g., genes, proteins, or other biological measurements).
* Columns represent the subjects or samples (e.g., individuals or experimental conditions).

For example, the rows of omics1 could represent genes, and rows of omics2 may be DNA copy number values.

Here’s how to replicate the steps above for real data:

- Simulate Omics Data: Imagine you’ve already preprocessed your omics data (omics1 and omics2).
- Get CT parameters: this step is identical to case of simulated data.
- Define feature set lists: featureSet1 and featureSet2 are lists of pathways of interest. 
- Estimate TDPs: this step is also identical to case of simulated data.

__More on feature set lists__: 
In omics studies, a feature set refers to a collection of genes, proteins, or other molecules that work together to carry out a specific biological function or process. These feature set represent interconnected biochemical reactions, molecular interactions, or regulatory networks within a cell, tissue, or organism. For example, in genomics, a feature set may refer to a set of genes involved in a particular process like cell cycle regulation, apoptosis (programmed cell death), or immune response. In proteomics, feature sets often involve proteins and their interactions in signaling cascades, metabolic cycles, or other molecular functions. And finally for DNA CN data, chromosome arms, cytobands, or known CNV hotspots are some common feture sets.\\
In the context of rOCEAN, we use predefined lists of features that belong to a particular feature set. These feature sets are often curated from databases like KEGG, Reactome, or MSigDB, which provide comprehensive collections of feature sets for various organisms. Each list should contain set of feature that correspond to a specific biological pathway. The eleents of the list are subsets of row names from the corresponding omics datasets. For example, featureSet1[[1]] might be a list of gene names associated with a specific Hallmark pathway. Combination of these feature sets will define the two-way feture sets. Refer to  [rSEA's vignette](https://github.com/cran/rSEA/blob/master/vignettes/rSEA_vignette.Rmd)] for more details on the feature set lists.

Now lets see some example code:

```{r}
#load omics datasets
load("PATH_TO_DATA/omics1.RData")
load("PATH_TO_DATA/omics2.RData")

#get closed testing parameters for the omics pair
gCT<-simesCT(omics1, omics2, alpha=0.05)
```
Once more, note that this is only done once for each pair of omics data.\\
Depending on the size of your data and your hardware, this step may take considerable time to complete. For large datasets, consider using a pre-computed $p$-value matrix and applying threshold.
If pvalmat is the pre-computed $p$-value matri, then using the following can speed up the calculation.

```{r}
spvalmat=pvalmat[pvalmat<0.05]
gCT<-simesCT(mps=pvalmat, m=nrow(pvalmat)*ncol(pvalmat)
```
Note that if the thresholded matrix is used, m is non-optional to pass the original dimention of the data to function.\\

The next step is estimation of TDPs.

```{r}
#load feature sets
load("PATH_TO_LISTS/featureSet1.RData")
load("PATH_TO_LISTS/featureSet2.RData")

#define the two-way feature set
feat_ids1<-which(rownames(omics1) %in% featureSet1[[1]])  # Pathway 1 in omics1
feat_ids2<-which(rownames(omics2) %in% featureSet2[[5]])  # Pathway 5 in omics2

#run ocean
out<-ocean(om1 = omics1[feat_ids1, ],
            om2 = omics2[feat_ids2, ],
            gCT = gCT, 
            scale = c("pair", "row", "col"))

```

This step can be repeated several times for all feature sets of interest.
