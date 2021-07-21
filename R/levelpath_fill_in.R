#' Find LevelPaths with holes
#' @param flowpaths a network obj
#' @return vector of dicconected levelpaths
#' @importFrom sf as_Spatial st_geometry_type st_as_sf st_line_merge
#' @importFrom rgeos gLineMerge
#' @importFrom dplyr mutate
#' @export

find_disconnected_lp = function(flowpaths){
  SPDF =  as_Spatial(fl)

  rownames(SPDF@data) <- sapply(slot(SPDF, "lines"), function(x) slot(x, "ID"))
  
  tmp <- rgeos::gLineMerge(SPDF, byid = TRUE, id = fl$levelpath)
  
  ids <- sapply(slot(tmp, "lines"), function(x) slot(x, "ID"))
  
  out = st_as_sf(tmp) %>% 
    mutate(levelpath = ids) %>% 
    st_line_merge() 
  
  runner = out$levelpath[st_geometry_type(out2) == "MULTILINESTRING"]
  
  unique(runner) 
}


#' Fill disconnected LevelPaths 
#' @param lpID levelpath ID to fill
#' @param network network
#' @importFrom dplyr filter slice_min slice_max mutate bind_rows
#' @importFrom sf st_intersects st_union st_buffer st_line_merge st_cast st_as_sf st_collection_extract st_filter
#' @importFrom sfnetworks as_sfnetwork activate to_spatial_shortest_paths edge_length
#' @importFrom tidygraph edge_is_multiple edge_is_loop convert
#' @importFrom lwgeom st_split
#' @return
#' @export

fill_level_path = function(lpID, network){
  # Base level path
  lp = filter(network, levelpath == lpID) 
  # Find current top and tail
  og_head  = slice_min(lp,hydroseq)
  og_tail  = slice_max(lp,hydroseq)
  
  h1 = find_node(x = og_head, "start")
  h2 = find_node(og_head, "end")
  head = ifelse(lengths(st_intersects(h1, lp)) > lengths(st_intersects(h2, lp)),
                h2, h1)
  
  t1 = find_node(og_tail, "start")
  t2 = find_node(og_tail, "end")
  tail = ifelse(lengths(st_intersects(t1, lp)) > lengths(st_intersects(t2, lp)),
                t2, t1)
  
  candidate = st_filter(network, st_union(st_buffer(lp, 5000)))

  t = st_line_merge(st_union(candidate)) %>% 
    st_cast("LINESTRING") %>% 
    st_as_sf()

  
  net = as_sfnetwork(t, directed = FALSE) %>%
    activate("edges") %>%
    filter(!edge_is_multiple()) %>%
    filter(!edge_is_loop()) %>% 
    convert(
      to_spatial_shortest_paths,
      from = head[[1]], to = tail[[1]],
      weights = edge_length()
    )

  new_lp = st_as_sf(activate(net, "edges")) 
  
  st_anti_filter = function(.x, .y, .predicate = st_intersects) {
    filter(.x, lengths(.predicate(.x, .y)) == 0)
  }
  
  c = build_node_net(candidate)
  
  breaks = st_collection_extract(lwgeom::st_split(candidate, c$node),"LINESTRING") 
   
  new_lp_full = st_filter(breaks, new_lp, .predicate = st_covered_by) %>% 
    mutate(levelpath = lpID)
  
  frags = st_anti_filter(breaks,  new_lp, .predicate = st_within) %>% 
    st_filter(new_lp, .predicate = st_touches) %>% 
    filter(ID %in% new_lp_full$ID) 
  
  corrections = list()
  
  if(nrow(frags) > 0){
    for(j in 1:nrow(frags)){
      cand = frags[j,]
      
      non_lp_connection = st_filter(c$network, cand,  .predicate = st_touches) %>% 
        filter(!ID %in% new_lp_full$ID)

      non_lp_connection$geometry = build_flow_line(frags$geometry[j], 
                                                   non_lp_connection$geometry)
      
      non_lp_connection$comid = paste(non_lp_connection$comid,
                                             frags$comid[j], sep = ",")
      
      
      corrections[[j]] = non_lp_connection
    } 
  } else {
    corrections = NULL
  }
  
  corr = bind_rows(corrections) 
  
  new_lp_full = new_lp_full %>% 
    mutate(levelpath = as.numeric(levelpath))
  
  improvements = bind_rows(new_lp_full, corr) 
  
  bind_rows(improvements, filter(network, !ID %in% improvements$ID))
}

