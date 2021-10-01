## ---------------------------
## Prepare, refactor, aggregate, and annotate a hydrofabric 
## for assigned RPU
##
## Author: Mike (mikecp11@gmail.com)
##
## ---------------------------
## Notes:
##
## ---------------------------

# Move down a level in project
setwd("workspace")
# Prettify the prints
options(paint_remind_mask_print = FALSE)

## ------------------------------------------------------

# Load creds and packages
library(rmarkdown)
library(paint); mask_print()
source("../inst/aws.R") # Personal file containing AWS access

# -------------------------------------------------------------------------
# User Defined Parameters
rpu_code   <- '01a'

data_dir   <- "NEEDED-DATA"
output_dir <- "OUTPUT"
AWS_bucket <- "formulations-dev/hydrofabric"

version    <- "beta"

meta <- read.csv('https://raw.githubusercontent.com/mikejohnson51/ngen-refactor-params/master/ngen-hydrofabric-mapping.csv')

cores      <-  2
# GO! ---------------------------------------------------------------------

(input <-  meta[meta$rpu_code == rpu_code, ])

## Step 1: NHD-navigate

render(
  input  = '01_NHD_navigate.Rmd',
  params = list(
    RPU        = input$rpu_code,
    min_da_km  = input$rf_min_da_km,
    data_dir   = data_dir,
    output_dir = output_dir,
    AWS_bucket = file.path(AWS_bucket, "01-reference-fabric")),
  envir        = new.env(),
  output_file  = paste0('temp/01_log_', rpu_code, '.html')
)

## Step 2: Refactor
### ASSUMES: reference fabric exists in output directory
render(
  input  = "02_refactor.Rmd",
  params = list(
    RPU                         = input$rpu_code,
    min_da_km                   = input$rf_min_da_km,
    data_dir                    = data_dir,
    output_dir                  = output_dir,
    AWS_bucket = file.path(AWS_bucket, "02-refactor-fabric"),
    #####
    split_flines_meters         = input$rf_split_flines_meters,
    cores                       = cores,
    collapse_flines_meters      = input$rf_collapse_flines_meters,
    collapse_flines_main_meters = input$rf_collapse_flines_main_meters
  ),
  envir       = new.env(),
  output_file = paste0('temp/02_log_', input$rpu_code, '.html')
)


# Step 3: Aggregate Catchments
render(
  input  = "03_aggregate.Rmd",
  params = list(
    RPU           = input$rpu_code,
    output_dir    = output_dir,
    AWS_bucket    = file.path(AWS_bucket, "03-aggregate-fabric"),
    
    keep          = input$agg_keep,
    ideal_area_km = input$agg_ideal_size_km,
    min_area_km   = input$agg_min_size_km,
    min_length_km = input$agg_min_length_km
  ),
  envir       = new.env(),
  output_file = paste0('temp/03_log_', input$rpu_code, '.html')
)

# Step 4: Release
render(
  input  = "04_release.Rmd",
  params = list(
    RPU      = input$rpu_code,
    dir      = release_dir,
    data_dir = output_dir,
    AWS_bucket = file.path(AWS_bucket, "releases"),
    version  = "beta"
  ),
  envir       = new.env(),
  output_file = paste0('temp/04_log_', input$rpu_code, '.html')
)
