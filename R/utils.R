flowpath_names = c("ID", "levelpath", "hydroseq", "comids",  "geometry")
catchment_names = c("ID", 'areasqkm', "geometry")

#' needs_layer
#' @description Checks if layer existing in geopackage
#' @param db character geopackage to check
#' @param layer character layer name
#' @return logical
#' @export
#' @importFrom sf st_layers


needs_layer <- function(db, layer) {
  
  if(file.exists(db)) {
    layers <- st_layers(db)
    if(layer %in% layers$name)
      return(FALSE)
  }
  TRUE
}

refactor_wrapper = function(flowpaths, 
                            catchments,
                            events = NULL,
                            avoid = NULL,
                            split_flines_meters = 10000, 
                            collapse_flines_meters = 1000,  
                            collapse_flines_main_meters = 1000,
                            cores = 1,  
                            facfdr = NULL,
                            routing = NULL,
                            keep = .9,
                            outfile){
  
  tf <- tempfile(pattern = "refactored", fileext = ".gpkg")
  tr <- tempfile(pattern = "reconciled", fileext = ".gpkg")
  
  if(!is.null(events)){
    events = filter(events, COMID %in% flowpaths$COMID)
  }
  
  if(!is.null(events)){
    avoid = avoid[avoid %in% flowpaths$COMID]
  }
  
  refactor_nhdplus(nhdplus_flines              = flowpaths, 
                   split_flines_meters         = split_flines_meters, 
                   split_flines_cores          = 1, 
                   collapse_flines_meters      = collapse_flines_meters,
                   collapse_flines_main_meters = collapse_flines_main_meters,
                   out_refactored = tf, 
                   out_reconciled = tr, 
                   three_pass          = TRUE, 
                   purge_non_dendritic = FALSE, 
                   events = events,
                   exclude_cats = avoid,
                   warn = FALSE)
  
  rec = st_transform(read_sf(tr), 5070)
  
  if(!is.null(routing)){
    
    rec$order = nhdplusTools::get_streamorder(st_drop_geometry(select(rec, ID, toID)), status = FALSE)
    
    rec = hyRefactor::add_lengthmap(rec, nhdplusTools::get_vaa("lengthkm")) %>% 
      attributes_for_flowpaths(
        weight_col = "lengthMap",
        length_weight = TRUE,
        rl_vars = c(
          "link", "Qi", "MusK", "MusX", "n",
          "So", "ChSlp", "BtmWdth",
          "time", "Kchan", "nCC",
          "TopWdthCC", "TopWdth", "alt"),
        rl_path  = routing
      )
    
    
    write_sf(st_transform(rec, 5070), outfile, "refactored_flowpaths", overwrite = TRUE)
    
  } else {
    
    write_sf(st_transform(rec, 5070), outfile, "refactored_flowpaths", overwrite = TRUE)
    
  }
  
  
  if(!is.null(facfdr)){
    
    rpus = unique(flowpaths$RPUID)
    rpus = rpus[!is.na(rpus)]
    
    fdrfac_files = list.files(facfdr, pattern = rpus, full.names = TRUE)
    fdr = raster::raster(grep("_fdr", fdrfac_files, value = TRUE))
    fac = raster::raster(grep("_fac", fdrfac_files, value = TRUE))
    catchments <-  st_transform(catchments, st_crs(fdr)) 
    st_precision(catchments) <- raster::res(fdr)[1]
    
    if("featureid" %in% names(catchments)){
      catchments = rename(catchments, FEATUREID = featureid)
    }
    
    reconciled <- st_transform(read_sf(tr), st_crs(fdr)) 
    refactored <- st_transform(read_sf(tf),  st_crs(fdr)) 
    
    divides    <- reconcile_catchment_divides(catchment = catchments,
                                              fline_ref = refactored,
                                              fline_rec = reconciled,
                                              fdr       = fdr,
                                              fac       = fac,
                                              para      = cores, 
                                              cache     = NULL, 
                                              fix_catchments = TRUE) 
    
    write_sf(st_transform(divides, 5070), outfile, "refactored_catchments", overwrite = TRUE)
  } 
  
  unlink(list(tr, tf))
  return(outfile)
}

attributes_for_flowpaths = function(flowpaths,
                                    weight_col = "lengthMap",
                                    length_weight = TRUE,
                                    rl_vars,
                                    rl_path){
  
  if(!"Length" %in% rl_vars){ rl_vars = c("Length", rl_vars) }
  
  net_map  <- dplyr::select(st_drop_geometry(flowpaths), ID, !!weight_col) %>%
    mutate(comid = strsplit(get(weight_col), ",")) %>%
    tidyr::unnest(cols = comid) %>%
    mutate(full_comids = floor(as.numeric(comid)),
           w = 10 * (as.numeric(comid) - full_comids), 
           w = ifelse(rep(length_weight, n()), w, 1),
           comid = NULL) 
  
  nc = RNetCDF::open.nc(rl_path)
  
  ll = lapply(rl_vars, function(x) RNetCDF::var.get.nc(nc, x))
  
  df = data.frame(do.call(cbind, ll)) %>% 
    setNames(rl_vars) %>% 
    rename(comid = link) %>% 
    right_join(net_map, by = c('comid' = 'full_comids')) %>% 
    select(-!!weight_col) %>% 
    mutate(w = w * Length) %>% 
    group_by(ID) %>% 
    summarise(across(everything(), ~ round(
      weighted.mean(x = ., 
                    w = w, 
                    na.rm = TRUE), 3))) %>% 
    dplyr::select(-comid, -Length, -w) 
  
  df2 = lapply(c("link", "gages", 'NHDWaterbodyComID'), function(x) x = RNetCDF::var.get.nc(nc, x)) %>% 
    bind_cols() %>% 
    setNames(c("link", "gages", 'NHDWaterbodyComID')) %>% 
    rename(comid = link) %>% 
    right_join(net_map, by = c('comid' = 'full_comids')) %>%
    mutate(gages = trimws(gages),
           gages = ifelse(gages == "", NA, gages),
           NHDWaterbodyComID = ifelse(NHDWaterbodyComID == -9999, NA, NHDWaterbodyComID)
    ) %>% 
    group_by(ID) %>% 
    summarise(gages = paste(gages[!is.na(gages)], collapse = ","),
              NHDWaterbodyComID = paste(unique(NHDWaterbodyComID[!is.na(NHDWaterbodyComID)]), collapse = ",")) %>% 
    left_join(df) %>% 
    mutate(gages = ifelse(gages == "", NA, gages),
           NHDWaterbodyComID = ifelse(NHDWaterbodyComID == "", NA, NHDWaterbodyComID))
  
  left_join(flowpaths, df2, by = "ID") %>% 
    mutate(Length_m = st_length(.))
}