#' Write GeoJSON
#' @param x sf object
#' @param file file path
#' @export
#' @importFrom sf write_sf st_make_valid st_transform
write_geojson <- function(x, file) {
  names(x) <- tolower(names(x))
  unlink(file)
  write_sf(st_make_valid(st_transform(x, 4326)), file, 
           layer_options = c("ID_FIELD=id", "ID_TYPE=String"))
}

build_toID_mapping = function(flowpaths){
  
  outlets <- flowpaths
  st_geometry(outlets) = st_geometry(find_node(outlets, "end"))
  
  imap = st_intersects(outlets, flowpaths)
  
  intmap = data.frame(ID = rep(outlets$ID, times = lengths(imap)),
                      toID = flowpaths$ID[unlist(imap)],
                      toHS = flowpaths$hydroseq[unlist(imap)]) %>% 
    filter(ID != toID) %>% 
    group_by(ID) %>% 
    slice_min(toHS, n = 1, with_ties = FALSE) %>% 
    ungroup() %>% 
    select(ID, toID)
  
  data.frame(
    ID = filter(outlets, !ID %in% intmap$ID)$ID,
    toID = 0) %>% 
    bind_rows(intmap) %>% 
    arrange(ID)
}


attributes_for_flowpaths = function(flowpaths,
                                    rl_vars,
                                    rl_path){
  
  if(!"Length" %in% rl_vars){
    rl_vars = c("Length", rl_vars)
  }
  
  net_map  <- select(st_drop_geometry(flowpaths), ID, comids) %>%
    mutate(comids = strsplit(comids, ",")) %>%
    tidyr::unnest(cols = comids) %>%
    mutate(comids = floor(as.numeric(comids))) %>%
    rename(comid = comids)
  
  nc = RNetCDF::open.nc(rl_path)
  
  ll = lapply(rl_vars, function(x) RNetCDF::var.get.nc(nc, x))
  
  df = data.frame(do.call(cbind, ll)) %>% 
    setNames(rl_vars) %>% 
    rename(comid = link) %>% 
    right_join(net_map, by = 'comid') %>% 
    group_by(ID) %>% 
    summarise(across(everything(), ~ round(
      weighted.mean(.x, 
                    .data$Length, 
                    na.rm = TRUE), 3))) %>% 
    select(-comid, -Length)
  
  df = lapply(c("link", "gages"), function(x) RNetCDF::var.get.nc(nc, x)) %>% 
    bind_cols() %>% 
    setNames(c("link", "gages")) %>% 
    rename(comid = link) %>% 
    right_join(net_map, by = 'comid') %>% 
    mutate(gages = ifelse(gages == "               ", NA, gsub('       ', "", gages))) %>% 
    group_by(ID) %>% 
    summarise(gages = paste(gages[!is.na(gages)], collapse = ",")) %>% 
    left_join(df) %>% 
    mutate(gages = ifelse(gages == "", NA, gages))
  
  left_join(flowpaths, df, by = "ID")
}

