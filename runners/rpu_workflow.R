system.time({
  rpu           = "01a"
  files         = build_rpu_parquet(rpu = rpu, overwrite = FALSE)
  orphans       = find_orphan_flowpaths(files) 
  cleaned       = merge_orphan_flowpaths(orphans)
  network_full  = prep_network_graph(cleaned, files$cat, ID_col = "ID")
  runner        = find_disconnected_lp(network_full$fl)
  console       = consolidate_levelpaths(network_full, runner)
  network_full2 = prep_network_graph(console, network_full$cat, ID_col = "ID")
  write_network_gpkg(network_full2, paste0("data/", rpu, "-test2.gpkg"))
})

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


