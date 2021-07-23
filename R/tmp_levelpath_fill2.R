#' Count flowpath linestrings
#' @param fp `sf` object
#' @return `sf` object
#' @export
#' @importFrom sf st_union st_line_merge st_cast st_as_sf

count_fps_ls = function(fp){
  fp %>% 
    st_union() %>% 
    st_line_merge() %>% 
    st_cast("LINESTRING") %>% 
    st_as_sf()
}

#' Anti filter
#' @param .x `sf` object
#' @param .y `sf` object
#' @param .predicate spatial predicate
#' @return `sf` object
#' @export
#' @importFrom dplyr filter

st_anti_filter = function(.x, .y, .predicate = st_intersects) {
  filter(.x, lengths(.predicate(.x, .y)) == 0)
}

#' Overlaps & Within
#' @param .x `sf` object
#' @param .y `sf` object
#' @return
#' @export
#' @importFrom sf st_within st_overlaps
st_overlaps_within = function(.x, .y) {
  int = lengths(st_within(.x, .y)) + lengths(st_overlaps(.x, .y))
  .x[int >= 1,]
}   

#lpID = runner[3]

############################################################################


#' Fill disconnected LevelPaths 
#' @param lpID levelpath ID to fill
#' @param network network
#' @importFrom dplyr filter slice_min slice_max mutate bind_rows pull slice
#' @importFrom sf st_intersects st_filter st_union st_buffer st_line_merge st_geometry st_as_sf st_collection_extract
#' @importFrom sfnetworks as_sfnetwork to_spatial_subdivision to_spatial_simple st_network_paths activate
#' @importFrom tidygraph convert `%E>%`
#' @importFrom lwgeom st_split
#' @importFrom foreach `%dopar%` foreach
#' @return
#' @export

fill_level_path2 = function(lpID, network){
  # Base level path network
  lp = filter(network, levelpath == lpID) 
  # Find current top and tail
  og_head  = slice_min(lp,hydroseq)
  og_tail  = slice_max(lp,hydroseq)
  
  # Find true head node
  h1 = find_node(x = og_head, "start")
  h2 = find_node(og_head, "end")
  head = ifelse(lengths(st_intersects(h1, lp)) > 
                lengths(st_intersects(h2, lp)),
                h2, h1)
  
  # find true tail node
  t1   = find_node(og_tail, "start")
  t2   = find_node(og_tail, "end")
  tail = ifelse(lengths(st_intersects(t1, lp)) > 
                lengths(st_intersects(t2, lp)),
                t2, t1)
  
  # Define candidate network
  e2 = st_filter(network, st_union(st_buffer(lp, 5000))) %>%
    flowpaths_to_linestrings()
  
  # Create undirected morphed graph
  net =  as_sfnetwork(e2, directed = FALSE) %>% 
    # Construct a subdivision of the network by subdividing edges at each 
    # interior point that is equal to any other interior or boundary point 
    # in the edges table. Interior points in this sense are those points 
    # that are included in their linestring geometry feature but are not
    # endpoints of it, while boundary points are the endpoints of the 
    # linestrings. The network is reconstructed after subdivision such
    # that edges are connected at the points of subdivision.
    convert(to_spatial_subdivision) %>% 
    # Remove loop edges and/or merges multiple edges into a single edge. 
    # Multiple edges are edges that have the same source and target nodes 
    # (in directed networks) or edges that are incident to the same nodes 
    # (in undirected networks). When merging them into a single edge, 
    # the geometry of the first edge is preserved. 
    convert(to_spatial_simple)
  
  # Find nodes of all connected linestring elements of levelpath
  nav =  find_node(count_fps_ls(lp), "both") 
  
  routes = foreach(i = 1:(nrow(nav)-1), 
                   .packages = c("sfnetworks", "sf")) %dopar% {
    st_network_paths(
      net, 
      from = st_geometry(nav)[i], 
      to   = st_geometry(nav)[i + 1]
    ) 
  }
  
  idx <- unique(unlist(pull(do.call("rbind", routes), edge_paths)))
  
  new_lp <- net %E>% 
    slice(idx) %E>% 
    st_as_sf() 
  
  ttt = count_fps_ls(new_lp)
  
  while(nrow(ttt) > 1){
    xx = as_sfnetwork(ttt) %>% 
      activate('edges') %>% 
      st_as_sf() %>% 
      st_filter(lp, .predicate = st_covers) 
    
    ttt = count_fps_ls(xx)
  }
  
  new_lp = st_overlaps_within(e2,ttt)
  
  c = build_node_net(new_lp)
  
  breaks = st_collection_extract(lwgeom::st_split(new_lp, c$node),
                                 "LINESTRING") %>% 
    mutate(breakID = 1:n())
  
  frags = st_anti_filter(breaks,  ttt, .predicate = st_within) %>% 
    st_filter(ttt, .predicate = st_touches) %>% 
    filter(ID %in% breaks$ID) 
  
  new_lp = filter(breaks, !breakID %in% frags$breakID)
  
  if (nrow(frags) > 0) {
    corrections = 
      foreach(i = 1:nrow(frags), .packages = c("sf")) %dopar% {
              frag.here = frags[j, ]
                            
              non_lp_connection = st_filter(network, frag.here,
                                            .predicate = st_touches) %>%
                              filter(!ID %in% breaks$ID)
                            
              if (nrow(non_lp_connection) != 0) {
                              st_geometry(non_lp_connection) =
                                build_flow_line(st_geometry(frag.here),
                                                st_geometry(non_lp_connection))
                              
                              non_lp_connection$comid = paste(non_lp_connection$comid,
                                                              frag.here$comid,
                                                              sep = ",")
                            
                              non_lp_connection
                            } else {
                              NULL
                            }
      }
  }
    
  #   for(j in 1:nrow(frags)){
  #     frag.here = frags[j,]
  #     
  #     non_lp_connection = st_filter(network, 
  #                                   frag.here,  
  #                                   .predicate = st_touches) %>% 
  #       filter(!ID %in% breaks$ID)
  #     
  #     
  #    if(nrow(non_lp_connection) != 0){
  #      st_geometry(non_lp_connection) = 
  #        build_flow_line(st_geometry(frag.here), 
  #                        st_geometry(non_lp_connection))
  #      
  #      non_lp_connection$comid = paste(non_lp_connection$comid,
  #                                      frag.here$comid, 
  #                                      sep = ",")
  #      
  #      
  #      corrections[[j]] = non_lp_connection
  #    } else {
  #      remove_ids[[j]] = frag.here$ID
  #    }
  #  } 
  # } else {
  #   corrections = NULL
  # }
  
  rids = frags$ID[which(unlist(lapply(corrections, is.null)))]
  
  mods =  bind_rows(mutate(new_lp, levelpath = as.numeric(lpID)), 
                   corrections)
  
  bind_rows(mods, filter(network, !ID %in% mods$ID)) %>% 
    filter(!ID %in% rids)
}
