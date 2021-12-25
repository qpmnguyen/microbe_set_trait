library(targets)
library(tarchetypes)
source("R/db_preprocess.R")

path_to_db <- function() {
    output <- file.path("data", "condensed_species_NCBI.txt")
    return(output)
}

list(
    tar_target(r_db, path_to_db(), format = "file"),
    tar_target(r_proc_db, process_db()),
    tar_target(r_tlist, get_trait_list(r_proc_db)),
    tar_target(r_sp_sets, create_sets(r_tlist, r_proc_db, agg = FALSE)),
    tar_target(r_gn_sets, create_sets(r_tlist, r_proc_db, agg = TRUE, threshold = 0.95)),
    tar_rds(r_sp_save, saveRDS(r_sp_sets, file = file.path("data", "sets_sp.rds"))),
    tar_rds(r_gn_save, saveRDS(r_gn_sets, file = file.path("data", "sets_gn.rds")))
)
