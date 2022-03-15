
<!-- README.md is generated from README.Rmd. Please edit that file -->
# HT-MMIOW

<!-- badges: start -->
<!-- badges: end -->
HT-MMIOW is a hypothesis test approach for microbiome mediation using inverse odds ratio. This approach uses singular value decomposition and UMAP to reduce the dimensions of microbiome data and quantifies the indirect effect of multiple mediators using inverse odds weighting. We employ a permutation test to test if an indirect effect exists.

## Installation

You can install the development version of HTMMIOW from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
library(devtools)
install_github("yukamoro/HTMMIOW")
```

And the development version from [GitHub](https://github.com/) with:

## Usage

Provide HTMMIOW() with exposure, mediator, and outcome data. The user may specify the number of umap components and permutations for the test.

``` r
# library(HTMMIOW)
# data(sim_data)
# 
# set.seed(10)
# htmmiow_output = Perform_HTMMIOW(exposure = sim_data$sim_exposure,
#                                  outcome = sim_data$sim_outcome,
#                                  mediator = sim_data$sim_microbiome,
#                                  umapcomponents = 2,
#                                  B = 2000)
# 
# htmmiow_output$pk_2
```

## References
