## The generalized gamma is a flexible distribution that outperforms alternatives when modelling catch

Code for an analysis evaluating the accuracy and predictive performance of the generalized gamma distribution (GGD) in population index standardization fishery survey data.
The analysis includes (1) a cross-simulation experiment and (2) spatiotemporal index standardization for scientific trawl survey data from the Gulf of Alaska and British Columbia using the R package sdmTMB.

## Citation:
Dunic, J.C., Conner, J., Anderson, S.C., and, J.T. Thorson. (2025). The generalized gamma is a flexible distribution that outperforms alternatives when modelling catch. ICES Journal of Marine Science. [https://doi.org/10.1093/icesjms/fsaf040](https://doi.org/10.1093/icesjms/fsaf040)

## Folders:

- `R`: code to run the analyses
- `data`: raw data; mostly kept locally due to size
- `data-outputs`: data that have been manipulated and cached; includes the cleaned data used in the multi-stock analysis
- `figures`: generated figures

## Analysis instructions:
To recreate the analyses, the R files are ordered based on dependencies. (i.e., generally `02-*` files cannot be run until `00-*` and `01-*` have been run). 

`00-prep-survey-data.R`: requires the raw data, but the outputs from this script are included in `data` and `data-outputs` 

## R Packages:

All R packages used should be available on CRAN with the following exceptions:

```
install.packages("pak")
pak::pkg_install(c(
  "seananderson/ggsidekick",
  "pbs-assess/gfplot",
  "pbs-assess/sdmTMB" # available on CRAN, but suggest using the latest development version:
))
```
