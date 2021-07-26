#' Build RPU Parquet file
#' Relies on have NHDPlus RDS files local
#' @param rpu RPU to extract
#' @param rds_dir RDS files
#' @param outdir Location to save data
#' @param overwrite Should exisitng files be overwritten?
#' @return list to files
#' @export
#' @importFrom dplyr filter select mutate `%>%` inner_join
#' @importFrom sf st_transform st_drop_geometry
#' @importFrom sfarrow st_write_parquet

build_rpu_parquet = function(rpu = "01a", rds_dir = "/Users/mjohnson/nhd_rds", 
                             outdir = "./data", overwrite = FALSE){
  
  RPUID <- COMID <- Hydroseq <- LevelPathI <- Shape <- 
    ID <- FEATUREID <- NULL
  
  outfile     = file.path(outdir, paste0(rpu,'_nhdplus_flowline_update.parquet'))
  outfile_cat = file.path(outdir, paste0(rpu,'_nhdplus_catchment.parquet'))
  
  if(any(!file.exists(outfile), !file.exists(outfile_cat), overwrite)){
  
    in_fl <- readRDS(file.path(rds_dir, 'nhdplus_flowline_update.rds')) %>%
      filter(RPUID == rpu) %>% 
      select(ID = COMID, 
             hydroseq = Hydroseq,  
             levelpath = LevelPathI, 
             geometry = Shape) %>%
      mutate(comid = ID) %>% 
      st_transform(5070) %>% 
      flowpaths_to_linestrings()
  
    outfile = file.path(outdir, paste0(rpu,'_nhdplus_flowline_update.parquet'))
    suppressWarnings({ sfarrow::st_write_parquet(obj=in_fl, dsn=outfile) })
    message("Finished ", paste0(rpu,'_nhdplus_flowline_update.parquet'))
  
    in_cat <- readRDS(file.path(rds_dir, 'nhdplus_catchment.rds')) %>%
      select(ID = FEATUREID, geometry = Shape) %>%
      inner_join(st_drop_geometry(in_fl), by = 'ID') %>% 
      st_transform(5070)
  
    suppressWarnings({ sfarrow::st_write_parquet(obj=in_cat, dsn=outfile_cat) })
    message("Finished ", paste0(rpu,'_nhdplus_flowline_update.parquet'))
  }
  
  return(list(cat = outfile_cat, fl = outfile))
}

