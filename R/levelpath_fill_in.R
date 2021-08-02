stubborn_flowpaths = function(x){
  bool = (st_geometry_type(sf::st_geometry(x)) == "MULTILINESTRING")
  multis = x[bool, ]
  geoms = list()
  for(i in 1:nrow(multis)){
    geoms[[i]] =  suppressWarnings({ 
      st_union_merge(st_cast(multis[i,], "LINESTRING")) 
    })
  }
  
  m = st_set_geometry(st_drop_geometry(multis), 
                      st_as_sfc(do.call(rbind, geoms), crs = st_crs(x)))
  
  do.call(rbind, list(m, x[!bool, ]))
}

st_union_merge = function(x){
  st_union(x) %>% 
    st_line_merge()
}

st_erase = function(x, y) st_difference(x, st_union(st_combine(y)))

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


#' rgeos LINE Merge
#' @description Wayyyy faster then either data.table, or sf based line merging
#' @param lines lines to merge
#' @param ID ID to merge over
#' @return 'sf object
#' @importFrom sf as_Spatial st_geometry_type st_as_sf st_line_merge
#' @importFrom rgeos gLineMerge
#' @importFrom dplyr mutate
#' @importFrom methods slot
#' @export

spLineMerge = function(lines, ID){
  
  SPDF =  as_Spatial(lines)
  
  rownames(SPDF@data) <- sapply(slot(SPDF, "lines"), function(x) slot(x, "ID"))
  
  tmp <- rgeos::gLineMerge(SPDF, byid = TRUE, id = lines[[ID]])
  
  ids <- as.numeric(sapply(slot(tmp, "lines"), function(x) slot(x, "ID")))
  
  st_as_sf(tmp) %>% 
    mutate("{ID}" := ids) %>% 
    flowpaths_to_linestrings() 
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
  
  out = spLineMerge(flowpaths, "levelpath")
  
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
#' @importFrom rlang `:=`
#' @export

fill_level_path = function(lpID, network){
  
  levelpath <- hydroseq <- to <- from <- . <-
    m <- weight <- ID <- .tidygraph_edge_index <- 
    breakID <- fromType <- toType <- NULL
  
  # Base level path
  lp = filter(network, levelpath == lpID) 
  # Find current top and tail
  og_head  = slice_min(lp,hydroseq)
  og_tail  = slice_max(lp,hydroseq)
  
  h1 = find_node(x = og_head, "start")
  h2 = find_node(og_head, "end")
  head = ifelse(lengths(st_intersects(h1, lp)) > 
                lengths(st_intersects(h2, lp)),
                h2, h1)
  
  t1 = find_node(og_tail, "start")
  t2 = find_node(og_tail, "end")
  tail = ifelse(lengths(st_intersects(t1, lp)) > 
                lengths(st_intersects(t2, lp)),
                t2, t1)
  
  candidate = st_filter(network, st_union(st_buffer(lp, 2500))) %>% 
    select(-to,-from) %>% 
    mutate(m = ifelse(levelpath == lpID, 0, 10000000)) %>% 
    #Adding weight to prioritize routing on levelpath IDs
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
    mutate(levelpath = as.numeric(lpID)) %>% 
    group_by(ID, hydroseq, levelpath, comid) %>% 
    summarize() %>% 
    ungroup()

  # Ensure the path includes define head and tail
  if(sum(sf::st_touches(og_head, new_lp_full) %>% lengths(),
         sf::st_touches(og_tail, new_lp_full) %>% lengths()) != 2) {
    stop("Only partial levelpath found!")
  }
  
  c = build_node_net(select(candidate, -m, -weight))
  
  breaks = st_collection_extract(st_split(c$network, c$node),"LINESTRING") 
  
  # Break lines and identify those that are NOT within the new lp, and that 
  # TOUCH the new lp, but are NOT the new lp
  frags  = st_anti_filter(breaks, new_lp_full, .predicate = st_within) %>%
    st_filter(new_lp_full, .predicate = st_touches) %>% 
    filter(ID %in% new_lp_full$ID) 

  return(list(lp = new_lp_full, f = frags))
  
  ## Option 1:
  # hw = filter(frags, fromType == "hw")
  # non_hw = filter(frags, fromType != "hw")
  # 
  # # Assign the frags the to ID of the from Nexus
  # to_merge = non_hw %>%
  #   select(from) %>%
  #   left_join(select(st_drop_geometry(breaks), -from), by = c('from' = 'to'))
  # 
  # o = filter(breaks, ID %in% to_merge$ID) %>% 
  #   bind_rows(to_merge) %>%
  #   select(-from, -to, -breakID, -fromType, -toType) %>% 
  #   arrange(ID)
  # 
  # if(nrow(o) > 0){
  #   o2 = st_drop_geometry(o) %>%
  #     filter(!duplicated(.)) 
  #   
  #   o3 = left_join(mutate(spLineMerge(o, 'ID'),
  #                         ID = as.integer(ID)), o2, by = "ID")
  # } else { 
  #   o3 = NULL
  # }
  # 
  # o3 %>% 
  #   bind_rows(new_lp_full, hw) %>% 
  #   select(hydroseq, levelpath, comid, ID)
  
  # Option 2
  # if(nrow(frags) > 0){
  #   for(j in 1:nrow(frags)){
  #       frag.here = frags[j, ]
  #       
  #       non_lp_connection = st_filter(candidate, frag.here,
  #                                     .predicate = st_touches) %>%
  #         filter(!ID %in% new_lp_full$ID)
  # 
  #       if(nrow(non_lp_connection) != 0){
  #         
  #         if(nrow(non_lp_connection) > 1){
  #           if(frag.here$levelpath %in% non_lp_connection$levelpath){
  #            non_lp_connection  = filter(non_lp_connection, levelpath == frag.here$levelpath)
  #           } else {
  #             non_lp_connection = slice_min(non_lp_connection, hydroseq)
  #           }
  #         }
  #   
  #         sf::st_geometry(non_lp_connection) =
  #           build_flow_line(sf::st_geometry(frag.here),
  #                           sf::st_geometry(non_lp_connection))
  #         
  #         non_lp_connection$comid = paste(non_lp_connection$comid,
  #                                         frag.here$comid,
  #                                         sep = ",")
  #         
  #         corrections[[j]] = non_lp_connection
  #       } else {
  #         corrections[[j]] = NULL
  #       }
  #   }
  # }

  # bind_rows(new_lp_full, corrections) %>% 
  #   select(-.tidygraph_edge_index) %>% 
  #   mapview() + tst
  
  #bind_rows(mods, filter(network, !ID %in% mods$ID)) 
}


#' Consolidate LevelPaths
#' @param network a network list
#' @param runner set of levelpaths to fill and correct
#' @importFrom parallel detectCores 
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach `%dopar%`
#' @importFrom dplyr filter bind_rows group_by summarize
#' @importFrom stats setNames
#' @export
#' 
consolidate_levelpaths = function(network, runner){
  
  i <- ID <- hydroseq <- levelpath <- comid <- NULL
  
  nCores <- parallel::detectCores() - 1
  doParallel::registerDoParallel(nCores)
  out = foreach(i = 1:length(runner)) %dopar% {
   fill_level_path(lpID  = runner[i], 
                   network = network$fl)
  }

  frags = lapply(out, '[[', 2) 
  frags = bind_rows(frags) %>% 
    mutate(FID = rep(1:length(frags), time = sapply(frags, nrow)))
  
  lps   = lapply(out, '[[', 1) 
  lps   = bind_rows(lps) %>% 
    flowpaths_to_linestrings() %>% 
    mutate(FID = rep(1:length(lps), time = sapply(lps, nrow))) %>% 
    select(hydroseq, levelpath, comid, ID, FID)
  
  others = network$fl %>% 
    filter(!ID %in% frags$ID) %>% 
    filter(!ID %in% lps$ID)
  
  ####
  
 dups = lps$ID[duplicated(lps$ID)]
 # RM Duplicates that were falsly merged into MUTLILINESTRINGS
 x = filter(lps, ID %in% dups) %>% 
     filter(st_geometry_type(.) == "LINESTRING")
 # Find Remaining Dups
 x2 = filter(x, ID %in% x$ID[duplicated(x$ID)]) 
 x  = filter(x, !ID %in% x2$ID)
 
# Find if
 x2$FF =  lapply(1:nrow(x2), function(x){
   st_contains_properly(st_union_merge(filter(lps, levelpath == x2$levelpath[x])),
                        x2) %>% 
     lengths()
 }) %>% unlist()
 
 lp_to_erase_from = filter(x2, FF == 0)
 collect = list()
 for(i in 1:nrow(lp_to_erase_from)){
   #erase y from x
   collect[[i]] = st_erase(filter(lps, levelpath == lp_to_erase_from$levelpath[i]), 
                                  lp_to_erase_from[i,])
 }
 
 new = bind_rows(collect)
 
 fl = filter(lps, !levelpath %in% unique(new$levelpath)) %>% 
   bind_rows(new) %>% 
   filter(!ID %in% x$ID) %>% 
   bind_rows(x, others) 

  if(nrow(network$cat) == nrow(fl)){
    return(fl)
  } else {
    stop("Unsuccessfull levelpath filling")
  }
}

