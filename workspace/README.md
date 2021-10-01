# ngen hydrofabric workspace

The markdown files in this directory can be used in development using RStudio or run programatically.They include:

The processing steps need to be executed across RPUs of NHDPlusV2 via a mapping file found at `../data/ngen-hydrofabric-mapping.csv`

1. `01_NHD-navigate.Rmd`
  - This file is adopted from the `gfv.20` project to ensure a common starting (or reference) fabric that resolves to a common set of POIs.
  
2. `02_hyRefactor.Rmd`
  - Implements a network flowline and catchment refactoring based on parameters defined in `../data/ngen-hydrofabric-mapping.csv`
  
3. `03_aggregate.Rmd`
  - Implements a network aggregation based on  parameters defined in `../data/ngen-hydrofabric-mapping.csv`

4. `04_release.Rmd`
  - Cuts a release of a given versioned hydrofabric for use in NGEN applications

5. `05_catchment_traits.Rmd` 
  - Produces catchment traits found in CAMELS for the released hydrofabric
