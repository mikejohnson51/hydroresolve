
rpu = "01a"
files   = build_rpu_parquet(rpu = rpu, overwrite = FALSE)
orphans = find_orphan_flowpaths(files) 
cleaned = merge_orphan_flowpaths(orphans)
network = prep_network_graph(cleaned, files$cat, ID_col = "ID")
write_network_gpkg(network, paste0("data/", rpu, "-test.gpkg"))
                   
runner = find_disconnected_lp(network$fl)
tmpNet = network$fl

system.time({
  for(i in 1:6){
    tmpNet = fill_level_path(lpID  = runner[i], network = tmpNet)
    table(st_geometry_type(tmpNet))
    message(i)
  }
})

devtools::document()
devtools::load_all()
rcompendium::add_dependencies()
rcompendium::add_dependencies_badge()
rcompendium::add_r_depend()


# rcompendium::add_github_actions_check(overwrite = TRUE)
# rcompendium::add_github_actions_check_badge()
# rcompendium::add_license("MIT", "Mike", "Johnson")
# rcompendium::add_license_badge()
# rcompendium::add_lifecycle_badge('experimental')
# rcompendium::add_repostatus_badge("concept")


