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

#' Find LevelPaths with holes
#' @param flowpaths a network obj
#' @return vector of dicconected levelpaths
#' @importFrom sf as_Spatial st_geometry_type st_as_sf st_line_merge
#' @importFrom rgeos gLineMerge
#' @importFrom dplyr mutate
#' @importFrom methods slot
#' @export

find_disconnected_lp = function(flowpaths){
  
  SPDF =  as_Spatial(flowpaths)

  rownames(SPDF@data) <- sapply(slot(SPDF, "lines"), function(x) slot(x, "ID"))
  
  tmp <- rgeos::gLineMerge(SPDF, byid = TRUE, id = flowpaths$levelpath)
  
  ids <- sapply(slot(tmp, "lines"), function(x) slot(x, "ID"))
  
  out = st_as_sf(tmp) %>% 
    mutate(levelpath = ids) %>% 
    st_line_merge() 
  
  runner = out$levelpath[st_geometry_type(out) == "MULTILINESTRING"]
  
  unique(runner) 
}

#' Fill disconnected LevelPaths 
#' @param lpID levelpath ID to fill
#' @param network network
#' @importFrom dplyr select filter slice_min slice_max mutate bind_rows
#' @importFrom sf st_intersects st_union st_buffer st_filter st_union st_buffer st_as_sf st_collection_extract st_touches st_within st_length st_geometry
#' @importFrom sfnetworks as_sfnetwork to_spatial_subdivision to_spatial_simple to_spatial_shortest_paths activate
#' @importFrom tidygraph convert
#' @importFrom lwgeom st_split
#' @export

fill_level_path = function(lpID, network){
  
  levelpath <- hydroseq <- to <- from <- . <-
    m <- weight <- ID <- .tidygraph_edge_index <- NULL
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
  
  candidate = st_filter(network, st_union(st_buffer(lp, 2500))) %>% 
    select(-to,-from) %>% 
    mutate(m = ifelse(levelpath == lpID, 0, 100000)) %>% 
    mutate(weight = as.numeric(st_length(.)) * m)
  
  net =  as_sfnetwork(candidate, directed = FALSE) %>% 
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
    convert(to_spatial_simple)  %>% 
    convert(
      to_spatial_shortest_paths,
      from = head[[1]], to = tail[[1]]
    )
  
  new_lp_full = st_as_sf(activate(net, "edges")) %>%
    mutate(levelpath = as.numeric(lpID), m = NULL, weight = NULL)

  if(sum(sf::st_touches(og_head, new_lp_full) %>% lengths(),
         sf::st_touches(og_tail, new_lp_full) %>% lengths()) != 2) {
    stop("Only partial levelpath found!")
  }
  
  candidate = select(candidate, -m, -weight)
  c = build_node_net(candidate)
  
  breaks = st_collection_extract(st_split(candidate, c$node),"LINESTRING") 
  
  frags  = st_anti_filter(breaks,  new_lp_full, .predicate = st_within) %>% 
    filter(ID %in% new_lp_full$ID) 

  corrections = list()
  
  if(nrow(frags) > 0){
    for(j in 1:nrow(frags)){
        frag.here = frags[j, ]
        
        non_lp_connection = st_filter(candidate, frag.here,
                                      .predicate = st_touches) %>%
          filter(!ID %in% new_lp_full$ID)

        if(nrow(non_lp_connection) != 0){
          
          if(nrow(non_lp_connection) > 1){
            if(frag.here$levelpath %in% non_lp_connection$levelpath){
             non_lp_connection = filter(non_lp_connection, levelpath == frag.here$levelpath)
            } else {
              non_lp_connection = slice_min(non_lp_connection, hydroseq)
            }
          }
    
          sf::st_geometry(non_lp_connection) =
            build_flow_line(sf::st_geometry(frag.here),
                            sf::st_geometry(non_lp_connection))
          
          non_lp_connection$comid = paste(non_lp_connection$comid,
                                          frag.here$comid,
                                          sep = ",")
          
          corrections[[j]] = non_lp_connection
        } else {
          corrections[[j]] = NULL
        }
    }
  }

  bind_rows(new_lp_full, corrections) %>% 
    select(-.tidygraph_edge_index)
  
  #bind_rows(mods, filter(network, !ID %in% mods$ID)) 
}

