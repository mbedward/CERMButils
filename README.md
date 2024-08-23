
# CERMButils

<!-- badges: start -->
<!-- badges: end -->

A grab-bag of R functions to work with bushfire data, developed for the Centre
for Environmental Risk Management of Bushfires, University of Wollongong.

## Installation

``` r
# install.packages("remotes")
remotes::install_github("mbedward/CERMButils")
```

## Example

Derive fire progression polygons from a data set of fire extent polygons.

``` r
library(sf)
library(CERMButils)

dat_extent_by_time <- st_read("fire_extent_data.shp")

dat_progressions <- make_progressions(dat_extent_by_time,
                                      time_cols = "DATETIME")

```
