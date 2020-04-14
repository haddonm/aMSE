
<!-- README.md is generated from README.Rmd. Please edit that file -->

# aMSE

<!-- badges: start -->

<!-- badges: end -->

This packages up a new Abalone Management Strategy Evaluation framework.
It has a novel structure that is based around the spatial ideas of
populations within spatial management units (SMUs), within a region:

  - region - highest geogrpahical level. Is simply the totality of the
    spatial management units.

  - SMU - spatial management units. In Tasmania these would currently be
    the classical statistical blocks.

  - population - literally a population. These are the active components
    within the simulation framework. The dynamics of the simulation are
    based around the populations, although, with positive larval
    dispersal (the default) there is some dependency of neighboiuring
    populations on each other.

This conceptual structure has implications for how the simulations are
conditioned during the management strategy testing. While it is still
the case tha tthe general properties of each SMU are used to set the
scene for each SMU’s component populations, the advent of sufficient GPS
logger data now allows a more formal definition of the productivity of
each population. Currently, the modelling assumes that each population
is linearly arranged around the coastline, and this will allow, where
available, specific catch/yield data to be used to define the bounds of
each population within an SMU.

An important change from previous designs is that larval dispersal is
implemented explicitly rather than implicitly being held constant. This
alters the dynamics so that analytical equilibrium methods no longer
work and we need to resort to iterative approaches to equilibrium if the
larval dispersal rate is \> 0.0.

Of course, adding such a component increases the number of options that
may need to be explored, but currently we envisage including a very low,
a middle range, and a relatively high level of larval dispersal, to
examin its implications.

## Installation

You can install the development version from
[GitHub](https://github.com/haddonm/aMSE) with:

``` r
# install.packages("devtools")
devtools::install_github("haddonm/aMSE")
devtools::install_github(haddonm/rutilsMH) # also useful for plotting, etc
```

Alternatively, you can generate a branch that you can work on by cloning
the repository, which, again, can be done very simply within RStudio.
Open the New Project option in the project dialog at the top right of
the RStudio screen and selection Version Control, then use
‘<https://github.com/haddonm/aMSE>’ in the top box, identify where you
want the new directory put, and press return.

It would be a good idea to read Hadley Wickham’s draft chpater on Git
and GitHub at <https://r-pkgs.org/index.html>.

## Example

This is a basic example which illustrates the generation of an initial
equilibrium. It uses built in data-sets but usually you would generate
regionC (constants) and regionD (the dynamic components) from the
conditioning data:

``` r
# basic example code
library(aMSE)
library(rutilsMH)

  data(region1)
  glb <- region1$globals
  data(testregC)
  data(testregD)
  regDe <- testequil(regC=testregC,regD=testregD,glob=glb)
#> [1] matureB OK
#> [1] exploitB OK
#> [1] recruitment OK
#> [1] spawning depletion OK
  str(regDe)
#> List of 10
#>  $ matureB : num [1:40, 1:6] 1037 1037 1037 1037 1037 ...
#>  $ exploitB: num [1:40, 1:6] 1168 1168 1168 1168 1168 ...
#>  $ catch   : num [1:40, 1:6] 0 0 0 0 0 0 0 0 0 0 ...
#>  $ harvestR: num [1:40, 1:6] 0 0 0 0 0 0 0 0 0 0 ...
#>  $ cpue    : num [1:40, 1:6] NA 417 417 417 417 ...
#>  $ recruit : num [1:40, 1:6] 629906 629906 629906 629906 629906 ...
#>  $ deplsB  : num [1:40, 1:6] 1 1 1 1 1 ...
#>  $ depleB  : num [1:40, 1:6] 1 1 1 1 1 ...
#>  $ catchN  : num [1:105, 1:40, 1:6] 0 0 0 0 0 0 0 0 0 0 ...
#>  $ Nt      : num [1:105, 1:40, 1:6] 6.30e+05 2.40e-06 4.41e-05 6.75e-04 8.57e-03 ...
```
