system.time({
  rpu = "01a"
  files        = build_rpu_parquet(rpu = rpu, overwrite = FALSE)
  orphans      = find_orphan_flowpaths(files) 
  cleaned      = merge_orphan_flowpaths(orphans)
  network_full = prep_network_graph(cleaned, files$cat, ID_col = "ID")
  runner       = find_disconnected_lp(network_full$fl)
})
#write_network_gpkg(network, paste0("data/", rpu, "-test.gpkg"))
                   
catch = list()
system.time({
  for(i in 1:length(runner)){
    catch[[i]] = fill_level_path(lpID  = runner[i], network = network_full$fl)
    message(i)
  }
})

dim(network_full$cat)
mods = bind_rows(catch) 
tmp = bind_rows(mods, filter(network_full$fl, !ID %in% mods$ID)) 
dim(tmp)
table(tmp$ID) %>% sort(decreasing  = TRUE)

duplicated(st_geometry(filter(tmp, ID == 1763)))
mapview(filter(tmp, ID == 1763)[3,], zcol = "from")

devtools::document()
devtools::load_all()
devtools::check()
rcompendium::add_dependencies()
rcompendium::add_dependencies_badge()
rcompendium::add_r_depend()


# rcompendium::add_github_actions_check(overwrite = TRUE)
# rcompendium::add_github_actions_check_badge()
# rcompendium::add_license("MIT", "Mike", "Johnson")
# rcompendium::add_license_badge()
# rcompendium::add_lifecycle_badge('experimental')
# rcompendium::add_repostatus_badge("concept")


