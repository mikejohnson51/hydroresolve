# Move down a level in project
setwd("workspace")

## ------------------------------------------------------

# Load creds and packages
library(rmarkdown)
library(purrr)
library(dplyr)
source("../inst/aws.R") 

base         <- '/Volumes/Transcend/ngen'
ngen_mapping <-  'https://raw.githubusercontent.com/mikejohnson51/ngen-refactor-params/master/ngen-hydrofabric-mapping.csv'
nwm_dir      <- '/Volumes/Transcend/nwmCONUS-v216'

merge_refs = function(dir = file.path(base, "usgs-reference"), pattern = "03"){
  
  out_file = file.path(dirname(refs[1]), paste0("reference_", pattern, ".gpkg"))
  
  if(!file.exists(out_file)){
    refs = list.files(dir, pattern = pattern, full.names = TRUE)
    
    fps = lapply(1:length(refs), function(x) { read_sf(refs[x], 'nhd_flowline') })
    
    POIs = lapply(1:length(refs), function(x) { 
      poi_layer = paste0("POIs_", gsub(".gpkg", "", gsub("reference_", "", basename(refs[x]))))
      read_sf(refs[x], poi_layer) })
    
    fps  = bind_rows(fps)
    pois = bind_rows(POIs)
    
    write_sf(fps,  out_file, 'nhd_flowline')
    write_sf(pois, out_file,  paste0("POIs_", pattern))
  }
}

merge_refs(dir = file.path(base, "usgs-reference"), pattern = "03")
merge_refs(dir = file.path(base, "usgs-reference"), pattern = "10")

read.csv(ngen_mapping) %>% 
  mutate(VPU = substr(rpu_code, 1,2)) %>% 
  filter(as.numeric(VPU) <= 18) %>% 
  split(.$rpu_code) %>%
  map(~rmarkdown::render(
    input  = '01_ngen_reference.Rmd',
    params = list(
      RPU              =  .$rpu_code,
      min_da_km        =  .$rf_min_da_km,
      base_dir         =  base
      reference_fabric =  file.path(base, "usgs-reference"),
      output_dir       =  file.path(base, "CONUS-hydrofabric/ngen-reference"),
      AWS_bucket = "formulations-dev/hydrofabric/ngen-reference"),
    envir        = new.env(),
    output_file  = paste0('temp/01_ngen_reference_', .$rpu_code, '.html')
  ))

##########

## IF YOU WANT TO:
#   add catchment divides provide a directory to fdrfac tifs else NULL
#   add routing variables, add a path to a RouteLink file else leave NULL
#   upload to AWS add a bucket path and source your credentials else NULL

read.csv(ngen_mapping) %>% 
    mutate(VPU = substr(rpu_code, 1,2)) %>% 
    filter(as.numeric(VPU) <= 18) %>% 
    slice(4:8) %>% 
    split(.$rpu_code) %>%
    map(~rmarkdown::render(
      input  = '02_refactor_grid.Rmd',
      params = list(
        RPU              =  .$rpu_code,
        routelink        =  NULL,
        reference_fabric =  file.path(base, "CONUS-hydrofabric/ngen-reference"),
        output_dir       =  file.path(base, "CONUS-hydrofabric/ngen-refactor"),
        split_flines_meters = .$rf_split_flines_meters,
        collapse_flines_main_meters = .$rf_collapse_flines_main_meters,
        collapse_flines_meters = .$rf_collapse_flines_meters,
        fdrfac =  file.path(base, "fdrfac"),
        AWS_bucket = "formulations-dev/hydrofabric/ngen-refactor"),
      envir        = new.env(),
      output_file  = paste0('temp/02_ngen_refactor_', .$rpu_code, '.html')
    ))
})
