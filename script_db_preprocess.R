library(targets)
library(tarchetypes)
library(tidyverse)
library(glue)
source("R/db_preprocess.R")

db_path <- function() {
    output <- file.path("data", "condensed_species_NCBI.txt")
    return(output)
}


trait_list <- tibble(
    traits = c("pathways", "carbon_substrates", "sporulation", 
               "gram_stain", "cell_shape", "range_tmp", "range_salinity", 
               "motility")
)


file <- tar_target(raw_db_file, db_load(), format = "file")
db <- tar_target(db_file, fread(raw_db_file))
create_db <- tar_map(unlist = FALSE, values = )


list(
    file,
    db
)



list(
    tar_target(r_db, path_to_db(), format = "file"),
    tar_target(r_proc_db, process_db()),
    tar_target(r_tlist, get_trait_list(r_proc_db)),
    tar_target(r_sp_sets, create_sets(r_tlist, r_proc_db, agg = FALSE)),
    tar_target(r_gn_sets, create_sets(r_tlist, r_proc_db, agg = TRUE, threshold = 0.95)),
    tar_rds(r_sp_save, saveRDS(r_sp_sets, file = file.path("data", "sets_sp.rds"))),
    tar_rds(r_gn_save, saveRDS(r_gn_sets, file = file.path("data", "sets_gn.rds")))
)
