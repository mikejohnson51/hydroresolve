#' Find Orphan Flowlines
#'
#' @param files a list of filepaths (see `build_rpu_parquet`)
#' @importFrom sfarrow st_read_parquet
#' @importFrom sfnetworks as_sfnetwork activate 
#' @importFrom sf st_as_sf st_drop_geometry
#' @importFrom dplyr select filter left_join mutate group_by ungroup bind_rows anti_join
#' @export

#files = build_rpu_parquet(rpu = "01a", overwrite = FALSE)

find_orphan_flowpaths = function(files){
  
  in_fl   = sfarrow::st_read_parquet(files$fl)
  in_cat  = sfarrow::st_read_parquet(files$cat)
  
  network = as_sfnetwork(in_fl) %>%
    activate('edges') %>%
    sf::st_as_sf()
  
  meta = select(st_drop_geometry(network), ID, hydroseq, toLP = levelpath, comid)
  
  id_map = select(st_drop_geometry(network), from, to, ID, levelpath) %>%
    # FIND FL without CATS
    filter(!ID %in% in_cat$ID) %>%
    # Identify fromID from fromNODE
    left_join(select(st_drop_geometry(network), fromID  = ID, tmpto = to),
              by = c("from" = "tmpto")) %>%
    # Identify toID from toNODE
    left_join(select(st_drop_geometry(network), toID    = ID, tmpfrom = from),
              by = c("to" = "tmpfrom")) %>%
    # If to/from ID are not in cats, then set to NA
    mutate(toID = ifelse(toID %in% in_cat$ID, toID, NA),
           fromID = ifelse(fromID %in% in_cat$ID, fromID, NA)) %>%
    # Group By ID to determine unique inflows and outflow
    group_by(ID) %>%
    mutate(fromCount = sum(!is.na(fromID)),
           toCount   = sum(!is.na(toID))) %>%
    # Define first pass merge into
    mutate(send_to = ifelse(toCount == 1, toID, fromID)) %>%
    # drop
    select(-from, -to) %>%
    #ungroup
    ungroup()
  
  id_map2 = id_map %>%
    inner_join(meta, by = c("send_to" = "ID")) %>%
    group_by(ID) %>%
    mutate(id_count = n()) %>%
    ungroup()
  
  purgeable = anti_join(id_map, meta, by = c("send_to" = "ID"))

# -------------------------------------------------------------------------
  
  # Find IDs that appear twice and prioritize
  dup_id_keep = filter(id_map2, id_count > 1) %>% 
    filter(!ID %in% purgeable) %>% 
    mutate(dDS = levelpath == toLP) %>% 
    group_by(ID) %>% 
    slice_max(dDS) %>% 
    slice_min(hydroseq) %>% 
    ungroup() %>% 
    filter(!duplicated(select(.,ID, send_to))) %>% 
    select(-dDS)

  # Find cases where multiple IDs are being sent to the same toID
  multis = filter(id_map2, id_count == 1) %>%
    group_by(send_to) %>%
    mutate(n = n())   %>%
    ungroup() %>%
    filter(n > 1)

  # -------------------------------------------------------------------------
  # TRY reversing to_send to selection
  reverse_paths = multis %>%
    mutate(send_to = ifelse(fromCount == 1, fromID, toID)) %>%
    group_by(send_to) %>%
    mutate(n = n())

  #IF FROMs are NA and going to the same toID, drop the shorter one
  multis_stubborn = filter(reverse_paths, n > 1) %>%
    filter(sum(is.na(fromID)) != n) %>%
    ungroup()
  
  multis_rm = filter(reverse_paths, n > 1) %>%
    filter(sum(is.na(fromID)) == n) %>%
    ungroup()

  multis_keep = filter(reverse_paths, n == 1)

# -------------------------------------------------------------------------

  to_merge = id_map2 %>% 
    filter(!ID %in% reverse_paths$ID ) %>%
    filter(!ID %in% dup_id_keep$ID) %>% 
    bind_rows(multis_keep, multis_stubborn, dup_id_keep) %>% 
    select(-n, -fromCount, -toCount, - id_count, -toLP)
  
  to_purge = c(multis_rm$ID, purgeable$ID)

  orphans = list(network = network,
              to_merge = to_merge, 
              to_purge = to_purge)
  
  return(orphans)
}

#' Merge Orphan Flowlines
#' @param orphans a list of orphan flowpaths (see `find_orphan_flowpaths`)
#' @importFrom dplyr mutate summarize filter group_by bind_rows rename
#' @importFrom sf st_as_sf
#' @export

merge_orphan_flowpaths = function(orphans){

  to_merge = orphans$to_merge %>% 
    mutate(geom   = NA, 
           fromID = NULL,
           toID   = NULL)
  
  partial_net = filter(orphans$network, ID %in% c(to_merge$ID, to_merge$send_to))

  for(i in 1:nrow(to_merge)){
    here = filter(partial_net, ID == to_merge$ID[i])
    to   = filter(partial_net, ID == to_merge$send_to[i])
    
    to_merge$geom[i] = build_flow_line(here$geometry,
                                       to$geometry)
    to_merge$comid[i] = paste(here$comid, to$comid, sep = ",")
  }

  to_merge2 = to_merge %>%
    group_by(send_to) %>%
    mutate(n = n()) %>%
    st_as_sf(crs = 5070)

  singles = filter(to_merge2, n == 1) %>% 
    ungroup()

  combs   = filter(to_merge2, n  > 1) %>%
    group_by(send_to) %>%
    summarize(ID = ID[1],
              levelpath = levelpath[1],
              send_to = send_to[1],
              comid     = paste0(comid, collapse = ","),
              levelpath = levelpath[1],
              hydroseq  = hydroseq[1]) %>%
    flowpaths_to_linestrings()
  
  cleaned = bind_rows(singles, combs) %>%
    mutate(ID = NULL, n = NULL) %>%
    rename(ID = send_to, geometry = geom)
  
  cleaned$comid = sapply(strsplit(cleaned$comid, ","),
                       function(x) paste(unique(x), collapse = ","))

  
  og_net = orphans$network %>%
    filter(!ID %in% to_merge$ID) %>%
    filter(!ID %in% to_merge$send_to) %>%
    filter(!ID %in% orphans$to_purge) %>%
    mutate(comid = as.character(comid))

  bind_rows(og_net, cleaned)
}

