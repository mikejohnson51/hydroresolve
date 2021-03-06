---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->



# hydroresolve

<!-- badges: start -->
[![Dependencies](https://img.shields.io/badge/dependencies-13/50-red?style=flat)](#)
[![R CMD Check](https://github.com/mikejohnson51/hydroresolve/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/mikejohnson51/hydroresolve/actions/workflows/R-CMD-check.yaml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://choosealicense.com/licenses/mit/)
[![LifeCycle](man/figures/lifecycle/lifecycle-experimental.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![Project Status: Concept](https://www.repostatus.org/badges/latest/concept.svg)](https://www.repostatus.org/#concept)
<!-- badges: end -->


`Hydroresolve` is a workflow package built around refactoring, aggregating, and paramterizing/attributing hydrofabrics for `ngen`.

## Mapping File

A primary input to the `hydroresolve` workflow is a mapping file that specifies the refactoring and aggregation parameters to apply to each unit of the NHD This file is stored in the `./data/` directory and is called `ngen-hydrofabric-mapping.csv`

| Process  | Variable                       | Description   |
|----------|--------------------------------|---------------|
| ID       | RPU                            | Raster Processing Unit from the NHD  |
| refactor | rf_min_da_km                   | This parameter represents the contributing area for a catchment.Anything under that size is treated as a non-dendritic area. e.g. there is no network logic applied.   |
| refactor | rf_split_flines_meters         | The maximum length flowpath desired in the output.  |
| refactor | rf_collapse_flines_meters      | The minimum length of inter-confluence flowpath desired in the output.  |
| refactor | rf_collapse_flines_main_meters | The minimum length of between-confluence flowpaths.  |
| aggregate| agg_keep                       | Topologically-aware geometry simplification. `keep` defines the proportion of points to retain (0-1; default 0.05)  |
| aggregate| agg_ideal_area_k2             | The ideal catchment size in network  |
| aggregate| agg_min_area_km2              | The minimum catchment size in network  |
| aggregate| agg_min_length_km             | The minimum flowpath length allowed in network   |


```r
meta = read.csv('workspace/inst/ngen-hydrofabric-mapping.csv')
paint(meta)
#> data.frame [73, 9] 
#> rpu_code                       chr 01a 02a 02b 03a 03b 03c
#> rf_min_da_km                   int 20 20 20 20 20 20
#> rf_split_flines_meters         int 10000 10000 10000 10000 ~
#> rf_collapse_flines_meters      int 1000 1000 1000 1000 1000~
#> rf_collapse_flines_main_meters int 1000 1000 1000 1000 1000~
#> agg_keep                       dbl 0.9 0.9 0.9 0.9 0.9 0.9
#> agg_ideal_size_km              int 10 10 10 10 10 10
#> agg_min_size_km                int 3 3 3 3 3 3
#> agg_min_length_km              dbl 0.6 0.6 0.6 0.6 0.6 0.6
```


```r
t(meta[1,])
#>                                1      
#> rpu_code                       "01a"  
#> rf_min_da_km                   "20"   
#> rf_split_flines_meters         "10000"
#> rf_collapse_flines_meters      "1000" 
#> rf_collapse_flines_main_meters "1000" 
#> agg_keep                       "0.9"  
#> agg_ideal_size_km              "10"   
#> agg_min_size_km                "3"    
#> agg_min_length_km              "0.6"
```

## Data Needs

| Name     | Size                           | Type/read   | Description | 
|----------|--------------------------------|---------------|
| ID       | RPU                            | Raster Processing Unit from the NHD  |
| ID       | RPU                            | Raster Processing Unit from the NHD  |
| ID       | RPU                            | Raster Processing Unit from the NHD  |

| ID       | RPU                            | Raster Processing Unit from the NHD  |
| ID       | RPU                            | Raster Processing Unit from the NHD  |
| ID       | RPU                            | Raster Processing Unit from the NHD  |
| ID       | RPU                            | Raster Processing Unit from the NHD  |

| ID       | RPU                            | Raster Processing Unit from the NHD  |
| ID       | RPU                            | Raster Processing Unit from the NHD  |

| ID       | RPU                            | Raster Processing Unit from the NHD  |
| ID       | RPU                            | Raster Processing Unit from the NHD  |
| ID       | RPU                            | Raster Processing Unit from the NHD  |
