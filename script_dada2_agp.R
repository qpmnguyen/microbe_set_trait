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
    return(file.path(filepath))
}



list(
    tar_files(raw_data, list.files(file.path(path_files(), "raw"), 
                                    pattern = "_001.fastq", full.names = TRUE)),
    tar_target(filt_path, filter_and_trim(raw_data, base_path = path_files(), pat = "_001.fastq"),
               pattern = map(raw_data)),
    tar_files(filt_data, filt_path),
    tar_target(l_errors, learn_errors(filt_data)),
    tar_target(r_dada, run_dada2(filt_data, l_errors), map(filt_data), iteration = "list"),
    tar_target(r_chimera, remove_chimera(r_dada)), 
    tar_target(a_taxonomy, assign_taxonomy(r_chimera)),
    tar_rds(save_results, saveRDS(a_taxonomy, file = here("data", "agp_dada2.rds")))
)
