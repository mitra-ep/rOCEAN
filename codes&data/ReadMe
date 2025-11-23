# Real data example from paper

## This repository includes:

* Preprocessed TCGA datasets
* Closed testing (CT) parameter calculation
* Feature-set preparation scripts
* Analysis scripts

## Steps to recreate the results

### Prepare the data
The data are downloaded from the TCGA Research Network and have been preprocessed based on the steps described in Section 6 of the paper.
The following are the preprocessed data in .RData format:

- Breast cancer expression data: data_exp_breast_P1 & data_exp_breast_P2
(must be combined before analysis; split due to GitHub size limitations)
- Breast cancer DNA CN data: data_CN_breast
- Colorectal cancer expression data: data_exp_colon
- Colorectal cancer DNA CN data: data_CN_colon

### Create the closed testing parameters

As mentioned in the paper, closed testing parameters are unique for each dataset and depend on the feature sets used.
It is straightforward to calculate them using the `simesCT` function in the `rOCEAN` R-package.
Both the code and output are available:

* `runSimes_BC.Rmd`: CT parameters for breast cancer data, results saved as `gCT.RData`
* `runSimes_VC.Rmd`: CT parameters for colorectal cancer data, results saved as `gCT_colon.RData`

### Prepare the feature-set
Feature sets are created using different R packages and then refined to include only the probes from the relevant datasets listed above.
The codes for these steps are included in the `path_prep.Rmd` file. Note that these pathways are regularly updated, so the feature-set files used in the paper are also available here for reproducibility.
Below is a list of these pathway files along with a short description:

- `chrom.table.rda` and `chrom_table_updated.RData`: chromosome list along with chromosome arms and relevant locations
- `pathHall.RData`: Hallmark pathways refined for breast cancer expression data
- `pathArm.RData`: chromosome arms refined for breast cancer CN data
- `pathCB.RData`: chromosome bands refined for breast cancer CN data
- `pathHallcolon.RData`: Hallmark pathways refined for colorectal cancer expression data
- `pathArmcolon.RData`: chromosome arms refined for colorectal cancer CN data
- `pathCBcolon.RData`: chromosome bands refined for colorectal cancer CN dat

### Run the codes for each analysis

The analysis was performed on a cluster computer. The R scripts are designed such that the code is looped over every
two-way feature set, using the loop identifier "r," in the code.
The lists of two-way feature sets are saved separately to make the loop traceable:

- Breast cancer hallmark pathways by chromosome arm: `bc_arm.R`, looping over `allarm.RData`
- Breast cancer hallmark pathways by chromosome band: `bc_band.R`, looping over `allband.RData`
- Colorectal cancer hallmark pathways by chromosome arm: `cc_arm.R`, looping over `allarmcolon.RData`
- Colorectal cancer hallmark pathways by chromosome band (selected items): `cc_band.R`, looping over `allbandcolon.RData`



