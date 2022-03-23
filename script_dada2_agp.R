library(targets)
library(tarchetypes)
library(here)

options(tidyverse.quiet = TRUE)
tar_option_set(packages = c("dada2", "tidyverse", "here", "rlist"))

source(here("R", "dada2_func.R"))
base_path <- file.path("rc", "lab", "H", 
                       "HoenA", "Lab", "QNguyen", "ResultsFiles", "data", 
                       "agp")
# this function returns some paths (redefining it for the purposes of global def) 
path_files <- function(){
    base_path <- file.path("rc", "lab", "H", 
                           "HoenA", "Lab", "QNguyen", "ResultsFiles", "data", 
                           "agp")
    if (Sys.info()['sysname'] == "Darwin"){
        filepath <- file.path("/Volumes", base_path)
        if (!file.exists(file.path(filepath, ".placeholder"))){
            stop("This folder does not exist, please mount folder")
        } 
    } else if (Sys.info()['sysname'] == "Linux"){
        filepath <- file.path("/dartfs-hpc", base_path)
    } else {
        stop("OS not supported")
    }
    return(filepath)
}



list(
    tar_target(raw_data, path_files()),
    tar_target(f_data, filter_and_trim(raw_data, limit = FALSE)),
    tar_target(l_errors, learn_errors(f_data)),
    tar_target(r_dada, run_dada2(l_errors)),
    tar_target(a_taxonomy, assign_taxonomy(r_dada)),
    tar_rds(save_results, saveRDS(a_taxonomy, file = here("data", "agp_dada2.rds")))
)
