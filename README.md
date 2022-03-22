
<!-- README.md is generated from README.Rmd. Please edit that file -->
# HT-MMIOW

<!-- badges: start -->
<!-- badges: end -->
HT-MMIOW is a hypothesis test approach for microbiome mediation that uses inverse odds ratio weighting. In short, this approach uses singular value decomposition and UMAP to reduce the dimensions of microbiome data and quantifies the indirect effect of multiple mediators using inverse odds weighting. We employ a permutation test to test if an indirect effect exists.

## Installation

You can install the development version of HTMMIOW from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
library(devtools)
install_github("yukamoro/HTMMIOW")
```

## Usage

Provide `Perform_HTMMIOW` with exposure, mediator, and outcome data. The user may specify the number of umap components and permutations for the test.

``` r
library(HTMMIOW)
data(sim_data)

set.seed(10)
htmmiow_output = Perform_HTMMIOW(exposure = sim_data$sim_exposure,
                                 outcome = sim_data$sim_outcome,
                                 mediator = sim_data$sim_microbiome,
                                 umapcomponents = 2,
                                 B = 5000)

htmmiow_output$pk_2
#> [1] 0
```

Here, the p-value is 0, so we conlude that the microbiome is a mediator in the causal path between exposure and outcome. Further information can be found in the vignette and corresponding manuscript.

## Issues

Please report any issues [here](https://github.com/yukamoro/HTMMIOW/issues).

## Reference

Manuscript for HT-MMIOW available soon.
