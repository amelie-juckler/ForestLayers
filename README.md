
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ForestLayers

<!-- badges: start -->

<!-- badges: end -->

Automatic Stratification of Forest Vertical Layers

## Overview

*ForestLayers* is an R package currently in development. It is an
automated tool for parameterizing and visualizing vegetation vertical
structure from either 1D or 3D data on vegetation density. No advanced
programming skills are required, and it is designed for scientists
looking to apply new analytic methods and foresters seeking to explore
their forestry data. A publication with the description of the method is
in progress, the package is subject to future modifications.

## Installation

You can install the latest version of *ForestLayers* as follows:

``` r
if (!require("remotes")) install.packages("remotes")
devtools::install_github("amelie-juckler/ForestLayers")
```

## Help section

You can look at the documentation directly from RStudio with those
requests:

``` r
?VerticalLayers_from_1Dprofile
?VerticalLayers_from_3Dgrid
?export_DVV
```

## Utilization

*ForestLayers* operates on vertical vegetation profiles (1D) or
voxelized representations of the point cloud (3D). The package contains
two main functions: VerticalLayers_from_1Dprofile and
VerticalLayers_from_3Dgrid, and one export function: export_DVV.

The choice of main function depends solely on the format of the input
data; use VerticalLayers_from_1Dprofile with a 1D vertical profile and
VerticalLayers_from_3Dgrid with a voxel grid or a 3D matrix. These
functions consist of a series of internal functions that: normalize the
3D grid if necessary, convert it into a 1D vertical profile,
automatically identify the vertical strata, and fit a Weibull curve to
the vertical distribution of each stratum. They also allow you to create
a graph to visualize the fitting results.

The export function allows users to save their results in a clear,
structured format outside the R environment.

### Script example

``` r
library(ForestLayers)
```

In construction …

### Sample results

In construction …

## Scientific methodology

A detailed description of the methodology will be provided after the
paper’s publication.

## Sample data

In construction …
