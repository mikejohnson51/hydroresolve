system.time({
  rpu           = "01a"
  files         = build_rpu_parquet(rpu = rpu, overwrite = FALSE)
  orphans       = find_orphan_flowpaths(files) 
  cleaned       = merge_orphan_flowpaths(orphans)
  network_full  = prep_network_graph(flowpaths  = cleaned, 
                                     catchments = files$cat,
                                     ID_col = "ID",
                                     directed = TRUE)
  
  runner        = find_disconnected_lp(flowpaths = network_full$fl)
  
  console       = consolidate_levelpaths(network = network_full, runner)
  
  console = select(console, -FID)
  
  tmp = spLineMerge(console, "levelpath")
  
  table(st_geometry_type(tmp))

  network_full2 = prep_network_graph(console, 
                                     network_full$cat, 
                                     ID_col = "ID",
                                     directed = TRUE)
  
  write_network_gpkg(network_full2, paste0("data/", rpu, "-test2.gpkg"))
})



# create_release(fl = network_full2$fl,
#                cat = network_full2$cat,
#                dir = "release",
#                ver = "tmp")


fl = sf::read_sf(paste0("data/", rpu, "-test2.gpkg"), "flowpaths")
        

devtools::document()
devtools::load_all()
devtools::check()
rcompendium::add_dependencies()
rcompendium::add_dependencies_badge()
rcompendium::add_r_depend()
knitr::knit("README.Rmd")

# rcompendium::add_github_actions_check(overwrite = TRUE)
# rcompendium::add_github_actions_check_badge()
# rcompendium::add_license("MIT", "Mike", "Johnson")
# rcompendium::add_license_badge()
# rcompendium::add_lifecycle_badge('experimental')
# rcompendium::add_repostatus_badge("concept")


