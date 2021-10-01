flowpath_names = c("ID", "levelpath", "hydroseq", "comids",  "geometry")
catchment_names = c("ID", 'areasqkm', "geometry")

st_union_merge = function(x){
    st_line_merge(st_union(x))
}

#' Convert MULITLINESTINGS to LINESTRINGS
#' @param flowpaths a flowpath `sf` object
#' @return a `sf` object
#' @export
#' @importFrom sf st_geometry_type st_geometry st_line_merge
#' @importFrom dplyr bind_rows

flowpaths_to_linestrings = function(flowpaths){
  bool = (st_geometry_type(sf::st_geometry(flowpaths)) == "MULTILINESTRING")
  multis = flowpaths[bool, ]
  if(nrow(multis) > 0){
    sf::st_geometry(multis) = st_line_merge(sf::st_geometry(multis))
  }
  singles = flowpaths[!bool, ]
  
  bind_rows(multis, singles)
}

#' Build LINESTRING from 2 geoms
#' @param geom1 (MULTILINESTRING) object 1
#' @param geom2 (MULTILINESTRING) object 2
#' @return a `sf` object
#' @export
#' @importFrom sf st_geometry_type st_union st_line_merge

build_flow_line = function(geom1, geom2){
  g = st_union(c(geom1, geom2))
  
  if(st_geometry_type(g) == "MULTILINESTRING"){
    g = st_line_merge(g)
  }
  
  g
}

#' Read Hydro Object
#' @param obj either a parquet file or sf object
#' @return an sf object
#' @importFrom sfarrow st_read_parquet
#' @importFrom methods is
read_hydro = function(obj){
  
  if(is.character(obj)){
    if(grepl('parquet', obj)){
      return(sfarrow::st_read_parquet(obj))
    }
  } else if(methods::is(obj, "sf")) { 
    return(obj)
  } else {
    message("obj must be an sf object or parquet file path")
  }
}

#' Find flowline node
#' Returns the start or end node from a flowline.
#' @param x `sf` flowline object
#' @param position Node to find: "start" or "end" or "both"
#' @export
#' @importFrom sf st_coordinates st_as_sf st_crs
#' @importFrom dplyr group_by filter row_number ungroup

find_node = function (x, position = "end") {
  
  X <- Y <- L1 <- L2 <- NULL
  
  tmp <- as.data.frame(st_coordinates(x))
  
  if ("L2" %in% names(tmp)) {
    tmp <- group_by(tmp, L2)
  } else {
    tmp <- group_by(tmp, L1)
  }
  
  if (position == "end") {
    tmp <- filter(tmp, row_number() == n())
  } else if (position == "start") {
    tmp <- filter(tmp, row_number() == 1)
  } else {
    tmp <- filter(tmp, row_number() %in% c(1, n()))
  }
  
  tmp <- select(ungroup(tmp), X, Y)
  st_as_sf(tmp, coords = c("X", "Y"), crs = st_crs(x))
}

strip_reference_fabric = function(rpu_code, bucket = NULL){
  
  path       = file.path("cache", paste0(rpu_code, ".gpkg"))
  flowpaths  = read_sf(path, "reconciled")
  catchments = read_sf(path, "divides")
  
  sf::write_sf(flowpaths, 
               paste0("../data/", rpu_code, "-reference.gpkg"), 
               "flowpaths")
  
  sf::write_sf(catchments, 
               paste0("../data/", rpu_code, "-reference.gpkg"), 
               "catchments")
  
  if(!is.null(bucket)){
    aws.s3::put_object(
      file   = paste0("../data/", rpu_code, "-reference.gpkg"), 
      object = paste0(rpu_code, "-reference.gpkg"), 
      bucket = bucket,
      multipart = TRUE
    )
  }
}

add_area = function(x){
  as.numeric(units::set_units(st_area(x), "km2"))
}  

add_length = function(x){
  as.numeric(units::set_units(st_length(x), "km"))
}  

spCatMerge = function(poly, ID){
  SPDF =  as_Spatial(poly)
  rownames(SPDF@data) <- sapply(slot(SPDF, "polygons"), function(x) slot(x, "ID"))
  tmp <- rgeos::gUnaryUnion(spgeom = SPDF, id = poly[[ID]],
                            checkValidity = 0)
  ids <- as.numeric(sapply(slot(tmp, "polygons"), function(x) slot(x, "ID")))
  
  suppressWarnings({
    st_as_sf(tmp) %>%
      mutate("{ID}" := ids) %>%
      mutate(areasqkm = add_area(.)) %>% 
      st_cast("POLYGON") %>% 
      st_make_valid()
  })
  
}

reorder_segments <- function(parts){
  
  joints = st_intersection(parts)
  joints = joints[st_is(joints,"POINT"),]
  jgraph = igraph::graph_from_edgelist(do.call(rbind, 
                                               joints$origins), 
                                       directed=FALSE)
  
  s = find_node(x, "start")
  e = find_node(x, "end")
  
  o = sfnetworks::st_network_paths(net, from = s, to = e) %>% 
    pull(node_paths)
  
  net = sfnetworks::as_sfnetwork(parts)  
  
  sps = igraph::shortest_paths(jgraph, ends[1], ends[2])
  
  path = sps$vpath[[1]]
  return(parts[path,])
}

connect_multilinestrings = function(x){
  # Break MLS into LS
  parts =  ms_explode(x) %>% 
    mutate(n = 1:n())
  # Find end node(s)
  if(nrow(parts) > 2){stop("Logic only applies to single-gap MLS")}
  nodes = find_node(parts, position = "both")
  # Pick off upper triangle of dist matrix
  dist = st_distance(nodes)[1:2, 3:4]
  # Pick closest nodes
  ind = which(dist == min(dist), arr.ind = TRUE)
  # Build new LS from those nodes
  int = nodes[c(ind[1], ind[2] + 2),] %>% 
    st_union() %>% 
    st_cast("LINESTRING") 
  # Merge with other parts
  st_as_sf(c(int,parts$geometry)) %>% 
    st_union_merge() 
}

cs_group <- function(x, threshold) {
  cumsum <- 0
  group <- 1
  result <- numeric()
  for (i in 1:length(x)) {
    cumsum <- cumsum + x[i]
    if (cumsum > threshold) {
      group <- group + 1
      cumsum <- x[i]
    }
    result = c(result, group)
  }
  return (result)
}

spatial_topology_correction = function(gpkg, keep = NULL){

  in_cat <- read_sf(gpkg, layer = 'divides') %>%
    dplyr::select(ID) %>%
    st_set_crs(5070) %>%
    mutate(areasqkm = add_area(.))
  
  in_fl = read_sf(gpkg, layer = "reconciled") %>%
    dplyr::select(ID, toID, levelpath = LevelPathID, hydroseq = Hydroseq, 
           comids = member_COMID,
           geometry = geom) %>% 
    st_transform(5070) %>% 
    mutate(lengthkm = add_length(.)) %>% 
    flowpaths_to_linestrings()
  
  #TODO: quick fix for orphan segments
  in_fl = filter(in_fl, ID %in% in_cat$ID)
  
  #############################################################################
  ### FILL GAP LPS
  #############################################################################
  
  gaps = filter(in_fl, st_geometry_type(in_fl) == "MULTILINESTRING")
  
  new_geom = lapply(1:nrow(gaps), function(x) {
    connect_multilinestrings(gaps[x,])
  })
  
  st_geometry(gaps) =  do.call('c', new_geom)
  
  in_fl = filter(in_fl, !ID %in% gaps$ID) %>% 
    bind_rows(gaps)
  
  if(sum(st_geometry_type(in_fl) == "MULTILINESTRING") != 0){
    stop("MLs remain in the network")
  }
  
  #############################################################################
  
  cat <- in_cat %>%
    rmapshaper::ms_explode() %>%
    filter(!duplicated(.)) %>%
    mutate(area = as.numeric(st_area(.)))
  
  ids <- filter(cat, duplicated(cat$ID))
  
  cat_no_problem <- filter(cat, !ID %in% ids$ID)
  
  challenges = filter(cat, ID %in% ids$ID) %>%
    mutate(tmpID = 1:n())
  
  base_cats = challenges %>%
    group_by(ID) %>%
    slice_max(area) %>%
    bind_rows(cat_no_problem)
  
  fragments = filter(challenges, !tmpID %in% base_cats$tmpID)
  
  message(prettyNum(nrow(fragments),big.mark=",",scientific=FALSE)," fragments to clean...")
  
  frags = fragments %>%
    rmapshaper::ms_dissolve() %>%
    rmapshaper::ms_explode() %>%
    mutate(area = as.numeric(st_area(.))) %>%
    sf::st_make_valid()
  
  message(prettyNum(nrow(frags),big.mark=",",scientific=FALSE), " Consolidated fragments...")
  
  ints = suppressWarnings({
    st_intersection(frags, st_make_valid(base_cats)) %>%
      st_collection_extract("LINESTRING") %>%
      mutate(l = st_length(.)) %>%
      group_by(rmapshaperid) %>%
      slice_max(l, with_ties = FALSE)
  }) 
  
  tj = right_join(frags, 
                  dplyr::select(st_drop_geometry(ints), ID, rmapshaperid), 
                  by = "rmapshaperid") %>%
    bind_rows(base_cats) %>%
    dplyr::select(-rmapshaperid, -areasqkm, -tmpID) %>% 
    group_by(ID) %>%
    mutate(n = n()) %>% 
    rename(geometry = geom) %>% 
    ungroup()
  
  cat = suppressWarnings({ 
    spCatMerge(filter(tj, n > 1) , 'ID') %>% 
      st_cast("POLYGON") %>% 
      bind_rows(dplyr::select(filter(tj, n == 1), ID)) %>% 
      mutate(tmpID = 1:n())
  }) 
  
  dups = cat$ID[duplicated(cat$ID)]
  
  if(length(dups) > 0){
    message("Cleaning up ", 
            length(dups), 
            " duplicated catchment fragments")
    
    for(i in 1:length(dups)){
      here = filter(cat, ID == dups[i])
      
      tmap = st_intersects(here, 
                           filter(in_fl, ID == dups[i]))
      
      dissolve = here[lengths(tmap) == 0, ]
      
      opt = suppressWarnings({
        st_intersection(dissolve, cat) %>% 
          st_collection_extract("LINESTRING") %>%
          mutate(l = st_length(.)) %>%
          slice_max(l, with_ties = FALSE)
      })
      
      ind = which(cat$tmpID == opt$tmpID.1)
      
      cat$geometry[ind] = st_union(dissolve$geometry, cat$geometry[ind])
      cat = filter(cat, tmpID != dissolve$tmpID)
    }
    
    cat = dplyr::select(cat, -tmpID)
  }
  
  if(!is.null(keep)){
    message("Simplifying catchment boundaries: keep = ", keep)
    cat = rmapshaper::ms_simplify(cat, keep = keep, keep_shapes = TRUE)
  }
  
  list(flowpaths  = in_fl, 
       catchments = cat)
}

merge_levelpath = function(network_list = out, ideal_size = 10){
  
  # Preprocessing .... ------------------------------------------------------
  cat = network_list$catchments %>% 
    mutate(areasqkm = add_area(.))
  
  fl  = network_list$flowpaths %>%
    mutate(lengthkm = add_length(.)) %>% 
    left_join(st_drop_geometry(cat), by = "ID")
  
  network_group = fl %>%
    group_by(levelpath) %>% 
    arrange(-hydroseq) %>%
    mutate(ind = cs_group(areasqkm, ideal_size)) %>%
    ungroup()   %>% 
    group_by(levelpath, ind) %>% 
    mutate(n = 1:n())
  
  network_group$g = group_indices(network_group)
  network_group   = ungroup(network_group)
  
  catchment_group = left_join(cat, 
                              select(st_drop_geometry(network_group), 
                                     ID, n, g))
  
  meta = network_group %>%
    st_drop_geometry() %>% 
    group_by(g) %>%
    summarize(
      levelpath = unique(levelpath),
      hydroseq  = min(hydroseq),
      comids    = paste(comids, collapse = ",")
    )
  
  tmp = left_join(spLineMerge(network_group, "g"), meta)
  
  gaps = filter(tmp, st_geometry_type(tmp) == "MULTILINESTRING")
  
  if(length(gaps) > 0){
   new_geom = lapply(1:nrow(gaps), function(x) {
    st_cast(gaps$geometry[x], "MULTIPOINT") %>% 
      st_combine() %>% 
      st_cast("LINESTRING")
  })
  
  st_geometry(gaps)  =  do.call('c', new_geom)
  } else {
    gaps = NULL
  }
  
  tmp = filter(tmp, !g %in% gaps$g) %>% 
    bind_rows(gaps) %>% 
    rename(ID = g)
  
  if(sum(st_geometry_type(tmp) == "MULTILINESTRING") != 0){
    stop("MLs remain in the network")
  }
  
  cat_tmp = suppressWarnings({
    spCatMerge(catchment_group, "g") %>% 
      st_cast("POLYGON")
  })
  
  dups = cat_tmp$g[duplicated(cat_tmp$g)]
  
  if(lengths(dups) > 0){
    
    cat_tmp$tmpID = 1:nrow(cat_tmp)
    message("Cleaning up ", length(dups), " duplicated catchment fragments")
    
    for(i in 1:length(dups)){
      here = filter(cat_tmp, g == dups[i])
      tmap = st_intersects(here, filter(tmp, ID == dups[i]))

      dissolve = here[lengths(tmap) == 0, ]
      
      opt = suppressWarnings({
        st_intersection(dissolve, cat_tmp) %>% 
          st_collection_extract("LINESTRING") %>%
          mutate(l = st_length(.)) %>%
          slice_max(l, with_ties = FALSE)
      })
        
      ind = which(cat_tmp$tmpID == opt$tmpID.1)
      
      cat_tmp$geometry[ind] = st_union(dissolve$geometry, cat_tmp$geometry[ind])
      cat_tmp = filter(cat_tmp, tmpID != dissolve$tmpID)
    }
    
    cat_tmp = dplyr::select(cat_tmp, -tmpID)
  }
  
  cat_tmp = rename(cat_tmp, ID = g)
  
  list(flowpaths  = tmp,
       catchments = cat_tmp)
}

#' rgeos LINE Merge
#' @description Wayyyy faster then either data.table, or sf based line merging
#' @param lines lines to merge
#' @param ID ID to merge over
#' @return an sf object
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

remove_islands = function(network_list){
  
  cat = network_list$catchments
  fl  = network_list$flowpaths
  
  inters = st_intersects(cat)
  
  interior_rings = data.frame(ID = cat$ID, n = lengths(inters)) %>%
    filter(n == 2)
  
  int_list = inters[which(cat$ID %in% interior_rings$ID)]
  
  all = cat$ID[unlist(int_list)]
  
  exteriors = all[!all %in% interior_rings$ID]
  
  message("Dissolving ", nrow(interior_rings), " islands...")
  
  outs = list()
  for(i in 1:nrow(interior_rings)){
    outs[[i]] = cat[int_list[[i]],] %>%
      mutate(area = st_area(.), id = 1) %>%
      arrange(-area) %>%
      group_by(id) %>%
      summarize(ID = ID[1])
  }
  
  new_exteriors = bind_rows(outs) %>%
    select(-id) %>%
    group_by(ID) %>%
    summarise()
  
  new_cat = cat %>%
    filter(!ID %in% all) %>%
    bind_rows(new_exteriors) %>%
    mutate(areasqkm = as.numeric(st_area(.)/1e6))
  
  new_fl = fl %>%
    filter(ID %in% new_cat$ID) %>%
    inner_join(st_drop_geometry(new_cat), by = 'ID')
  
  list(flowpaths = select(new_fl, flowpath_names), 
       catchments = select(new_cat, catchment_names))
}




mainstem_reduce = function(network_list, 
                           min_area, 
                           min_length) {
  
  cat = network_list$catchments %>% 
    mutate(areasqkm = add_area(.))
  
  fl  = network_list$flowpaths %>%
    mutate(lengthkm = add_length(.),
           areasqkm = NULL) %>% 
    left_join(st_drop_geometry(cat), by = "ID")
  
  network  = fl %>% 
    st_drop_geometry() %>% 
    group_by(levelpath) %>% 
    mutate(toID   = lead(ID),
           fromID = lag(ID),
           flag = lengthkm < min_length | areasqkm < min_area) %>% 
    ungroup() 
  
  mid_reach = network %>% 
    filter(flag & !is.na(toID)) %>% 
    select(ID, into = toID) 
  
  outlet = filter(network, is.na(toID) & flag) %>% 
    select(ID, into = fromID) %>% 
    filter(!is.na(into))
  
  gg = bind_rows(mid_reach, outlet) %>% 
    group_by(into)
  
  gg$g = group_indices(gg)
  
  gg = gg %>% 
    ungroup() %>% 
    tidyr::pivot_longer(-g) %>% 
    select(ID = value,
           direction = name, 
           g) 
  
  tmp = left_join(gg, fl, by = 'ID') %>% 
    st_as_sf()
  
  meta = tmp %>%
    st_drop_geometry() %>% 
    group_by(g) %>%
    arrange(direction) %>% 
    summarize(
      ID = last(ID),
      levelpath = unique(levelpath),
      hydroseq  = last(hydroseq),
      comids    = paste(comids, collapse = ",")
    ) %>% 
    ungroup()
  
  tmp = left_join(spLineMerge(tmp, "g"), meta)  %>% 
    select(-g)
  
  new_fl = fl %>% 
    filter(!ID %in% gg$ID) %>% 
    bind_rows(tmp)
  
  gaps = filter(new_fl, st_geometry_type(new_fl) == "MULTILINESTRING")
  
  new_geom = lapply(1:nrow(gaps), function(x) {
    st_cast(gaps$geometry[x], "MULTIPOINT") %>% 
      st_combine() %>% 
      st_cast("LINESTRING")
  })
  
  st_geometry(gaps) =  do.call('c', new_geom)
  
  new_fl = filter(new_fl, !ID %in% gaps$ID) %>% 
    bind_rows(gaps)
  
  tmp_cat = left_join(gg, cat, by = 'ID') %>% 
    st_as_sf()
  
  tmp_cat = left_join(spCatMerge(tmp_cat, "g"), select(meta, g, ID))  %>% 
    select(-g)
  
  new_cat = cat %>% 
    filter(!ID %in% gg$ID) %>% 
    bind_rows(tmp_cat)
  
  dups = new_cat$ID[duplicated(new_cat$ID)]
  
  if(lengths(dups) > 0){
    
    new_cat$tmpID = 1:nrow(new_cat)
    message("Cleaning up ", length(dups), " duplicated catchment fragments")
    
    for(i in 1:length(dups)){
      here = filter(new_cat, ID == dups[i])
      tmap = st_intersects(here, filter(new_fl, ID == dups[i]))

      dissolve = here[lengths(tmap) == 0, ]
      
      if(nrow(dissolve)  == 0){
        here$areasqkm = add_area(here)
        dissolve = slice_min(here, areasqkm, n = 1)
      }
      
      opt = suppressWarnings({
        st_intersection(dissolve, new_cat) %>% 
          st_collection_extract("LINESTRING") %>%
          mutate(l = st_length(.)) %>%
          slice_max(l, with_ties = FALSE)
      })
      
      ind = which(new_cat$tmpID == opt$tmpID.1)
      
      new_cat$geometry[ind] = st_union(dissolve$geometry, new_cat$geometry[ind])
      new_cat = filter(new_cat, tmpID != dissolve$tmpID)
    }
    
  }
  
  list(flowpaths  = select(new_fl, flowpath_names), 
       catchments = select(new_cat, catchment_names))

}

remove_single_levelpaths = function(network_list,
                                    min_area   = 3, 
                                    min_length = .6){
  
  
  cat = network_list$catchments %>% 
    mutate(areasqkm = add_area(.))
  
  fl  = network_list$flowpaths %>%
    mutate(lengthkm = add_length(.),
           areasqkm = NULL) %>% 
    left_join(st_drop_geometry(cat), by = "ID")
  
  net_nodes = fl %>% 
    group_by(levelpath) %>% 
    mutate(n = n()) %>% 
    ungroup() %>% 
    filter(n == 1) %>% 
    filter(areasqkm < 3 | lengthkm < .6)
  
  st_geometry(net_nodes) = st_geometry(find_node(net_nodes, 'end'))
  
  imap = st_intersects(net_nodes, fl)
  
  int_map = data.frame(
    nodeID = rep(net_nodes$ID, times = lengths(imap)),
    nodeLENGTH = rep(net_nodes$lengthkm, times = lengths(imap)),
    nodeAREA   = rep(net_nodes$areasqkm, times = lengths(imap)),
    intID  = fl$ID[unlist(imap)],
    intHS  = fl$hydroseq[unlist(imap)]
  )  %>% 
    filter(nodeID != intID) %>% 
    group_by(nodeID) %>% 
    slice_min(intHS) %>% 
    ungroup() %>% 
    filter() %>% 
    filter(!nodeID %in% intID) %>% 
    select(nodeID, intID) %>% 
    filter(!duplicated(.)) %>% 
    group_by(intID)
  
  int_map$g = group_indices(int_map)
  
  int_map_p = ungroup(int_map) %>% 
    tidyr::pivot_longer(-g, values_to = 'ID') %>% 
    filter(!duplicated(.))
  
  meta = filter(int_map_p, name == "intID") %>% 
    filter(!duplicated(.))
  
  tmp_cat = left_join(int_map_p, cat, by = 'ID') %>% 
    st_as_sf()
  
  c  = spCatMerge(tmp_cat, "g")
  cc = left_join(c, meta) 
  
  dups = filter(int_map_p, ID %in% cc$ID[which(duplicated(cc$ID))])$g
  
  return_to_unmerged = filter(int_map, g %in% dups)
  return_to_unmerged = unique(c(return_to_unmerged$nodeID, 
                                return_to_unmerged$intID))
  
  cc2 = filter(cc, !g %in% dups) %>% 
    select(-g, -name)
  
  mapping = filter(int_map, !g %in% dups)
  
  cat2 = cat %>% 
    filter(!ID %in% mapping$nodeID) %>% 
    filter(!ID %in% cc2$ID) %>% 
    bind_rows(cc2)
  
  fl2 = filter(fl, !ID %in% mapping$nodeID)
  
  list(flowpaths = fl2, 
       catchments = cat2)
}
