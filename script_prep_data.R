library(targets)
library(tarchetypes)
source("R/gevers_preprocess.R")
source("R/hmp_preprocess.R")
source("R/utils.R")


gevers_path <- function(){
  return("data/gevers_dada2.rds")
}





list(
  tar_target(r_gevers_path, gevers_path(), format = "file"),
  tar_target(r_gevers_physeq, load_and_format(file_path = r_gevers_path)),
  tar_target(r_gevers_sp,{
    physeq <- tax_glom(r_gevers_physeq, taxrank = "Species")
    physeq <- reformat_species(physeq)
    return(physeq)
  }), 
  tar_target(r_gevers_gn, {
    physeq <- tax_glom(r_gevers_physeq, taxrank = "Genus")
    return(physeq)
  })
  
)