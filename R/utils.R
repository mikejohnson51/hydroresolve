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

#' Build HydroNetwork Graph
#' @param flowpaths a filepath to a `parquet` file or `sf` object
#' @param catchments a filepath to a `parquet` file or `sf` object
#' @param ID_col the desired ID column
#' @return
#' @export
#' @importFrom sfnetworks as_sfnetwork activate
#' @importFrom sf st_as_sf
#' @importFrom dplyr mutate rename

prep_network_graph = function(flowpaths, catchments, ID_col = "ID"){
  
  flowpaths  = read_hydro(flowpaths)
  catchments = read_hydro(catchments)
  
  if(nrow(flowpaths) != nrow(catchments)){
    stop("Must have one flowpath per catchment!!")
  }
  
  id_map = data.frame(
    old_id = sort(flowpaths[[ID_col]]),
    new_id = 1:nrow(flowpaths)
  )
  
  
  fl_net = suppressWarnings({ as_sfnetwork(flowpaths) })
  
  nodes = st_as_sf(mutate(activate(fl_net,"nodes"), nexID = 1:n()))
  
  fl    = st_as_sf(mutate(activate(fl_net,"edges"))) %>% 
    rename(old_id = !!ID_col) %>% 
    mutate(ID = id_map$new_id[match(old_id, id_map$old_id)],
           old_id = NULL)
  
  cat    = catchments  %>% 
    rename(old_id = !!ID_col) %>% 
    mutate(ID = id_map$new_id[match(old_id, id_map$old_id)],
           old_id = NULL)
  
  hw   = fl$from[!fl$from %in% fl$to]
  term = fl$to[!fl$to %in% fl$from]
  
  nodes$type = ifelse(nodes$nexID %in% hw, "hw", NA)
  nodes$type = ifelse(nodes$nexID %in% term, "term", nodes$type)
  nodes$type = ifelse(is.na(nodes$type), "nex", nodes$type)
  
  return(list(nex = nodes, fl = fl, cat  = cat))
}


#' Write a HydroNetwork GPKG
#'
#' @param network a newtwork list containing a cat, fl, and nex item (see `prep_network_graph`)
#' @param outpath where to write the file?
#' @return
#' @export
#' @importFrom sf write_sf

write_network_gpkg = function(network, outpath){

  if(sum(names(network) %in% c('nex', 'fl', 'cat')) != 3){
    stop("Must have nex, fl, and cat list elements in network")
  }
  
  unlink(outpath)
  write_sf(network$fl,  outpath, layer = "flowpaths")
  write_sf(network$cat, outpath, layer = "catchments")
  write_sf(network$nex, outpath, layer = "nexi")
  
}


#' Build Node Net
#' @param flowpaths a filepath to a `parquet` file or `sf` object
#' @param add_type 
#' @importFrom sfnetworks as_sfnetwork activate
#' @importFrom tidygraph edge_is_multiple edge_is_loop
#' @importFrom sf st_as_sf st_drop_geometry
#' @importFrom dplyr mutate n left_join select
#' @return list
#' @export

build_node_net = function(flowpaths, add_type = TRUE){
  
  flowpaths = read_hydro(flowpaths)
  
  net = suppressWarnings({ sfnetworks::as_sfnetwork(flowpaths)  })
  # network    = st_as_sf(mutate(activate(net,"edges"))) 
  network = net %>% 
    activate('edges') %>% 
    filter(!edge_is_multiple()) %>%
    filter(!edge_is_loop()) %>% 
    st_as_sf()
  
  nodes = st_as_sf(mutate(activate(net,"nodes"), nexID = 1:n()))
  
  if(add_type){
    id_map = network %>% 
      # Identify fromID from fromNODE
      left_join(select(st_drop_geometry(network), fromID  = ID, tmpto = to),   
                by = c("from" = "tmpto")) %>% 
      # Identify toID from toNODE
      left_join(select(st_drop_geometry(network), toID    = ID, tmpfrom = from),
                by = c("to" = "tmpfrom"))
    
    hw   = id_map$from[!id_map$from %in% id_map$to]
    term = id_map$to[!id_map$to %in% id_map$from]
    nodes$type = ifelse(nodes$nexID %in% hw, "hw", NA)
    nodes$type = ifelse(nodes$nexID %in% term, "term", nodes$type)
    nodes$type = ifelse(is.na(nodes$type), "nex", nodes$type)
    
    network$fromType = left_join(select(network, from), 
                                 st_drop_geometry(nodes), by = c('from' = "nexID"))$type
    
    network$toType = left_join(select(network, to), 
                               st_drop_geometry(nodes), by = c('to' = "nexID"))$type
  }
  
  list(node = nodes, network = network)
}

#' Find flowline node
#' Returns the start or end node from a flowline.
#' @param x `sf` flowline object
#' @param position Node to find: "start" or "end"
#' @return
#' @export
#' @importFrom sf st_coordinates st_as_sf st_crs
#' @importFrom dplyr group_by filter row_number ungroup

find_node = function (x, position = "end") {
  tmp <- st_coordinates(x) %>% as.data.frame()
  if ("L2" %in% names(x)) {
    tmp <- group_by(tmp, .data$L2)
  } else {
    tmp <- group_by(tmp, .data$L1)
  }
  
  if (position == "end") {
    tmp <- filter(tmp, row_number() == n())
  } else if (position == "start") {
    tmp <- filter(tmp, row_number() == 1)
  }
  
  tmp <- select(ungroup(tmp), X, Y)
  st_as_sf(tmp, coords = c("X", "Y"), crs = st_crs(x))
}
