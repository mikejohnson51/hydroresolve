data_paths <- jsonlite::read_json(file.path("cache", "data_paths.json"))

all_rpus <- sort(gsub("rpu_", "", names(data_paths$fdr)))

build_mapping_csv = data.frame(
  #ID
  rpu_code = all_rpus,
  
  #Refactoring Parameters
  rf_min_da_km                   = 20,
  rf_split_flines_meters         = 10000,
  rf_collapse_flines_meters      = 1000,
  rf_collapse_flines_main_meters = 1000,
  
  #Aggregation Parameters
  agg_keep          = 0.9,
  agg_ideal_size_km = 10,
  agg_min_size_km   = 3,
  agg_min_length_km = 0.6
)

data.table::fwrite(build_mapping_csv, 
          'inst/ngen-hydrofabric-mapping.csv',
          row.names = FALSE
          )
