
<!-- README.md is generated from README.Rmd. Please edit that file -->

# sasfunclust

<!-- badges: start -->
<!-- badges: end -->

The package sasfunclust implements the the sparse and smooth functional
clustering (SaS-Funclust) method proposed by Centofanti et al. 2021.
SaS-Funclust is a new method for clustering functional data that aims to
classify a sample of curves into homogeneous groups while jointly
detecting the most informative portions of domain. The mothod relies on
a general functional Gaussian mixture model whose parameters are
estimated by maximizing a log-likelihood function penalized with the
functional adaptive pairwise fusion penalty and a roughness penalty.
\#\# Installation

You can install the released version of sasfunclust from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("sasfunclust")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("fabiocentofanti/sasfunclust")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(sasfunclust)
#> 
#> Attaching package: 'sasfunclust'
#> The following object is masked from 'package:graphics':
#> 
#>     plot
#> The following object is masked from 'package:base':
#> 
#>     plot
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date.

You can also embed plots, for example:

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub!
