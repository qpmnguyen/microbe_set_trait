library(targets)
library(tarchetypes)
library(tidyverse)
library(glue)
library(future.callr)
library(future.batchtools)
library(future)


plan(batchtools_slurm, template = "batchtools.slurm.tmpl")

source("R/db_preprocess.R")

get_db_path <- function() {
    output <- file.path("data", "madin_proc.rds")
    return(output)
}


trait_list <- cross_df(list(
    traits = c("pathways", "carbon_substrates", "sporulation", 
               "gram_stain", "cell_shape", "range_tmp", "range_salinity", 
               "motility"),
    g_agg = c(TRUE, FALSE)
))

print(getDTthreads())

get_path <- tar_target(db_path, get_db_path(), format = "file")
load_db <- tar_target(db, readRDS(raw_db_file))
create_db <- tar_map(unlist = FALSE, values = trait_list, 
    tar_target(t_list, get_traits(trait_db = db, trait = traits)),
    tar_target(t_set, get_sets(ncbiid_list = t_list, 
                               trait_db = db, g_agg = g_agg, 
                               trait = trait, prop_thresh = 0.95), 
               resources = tar_resources(
                   future = tar_resources_future(
                       plan = tweak(
                           batchtools_slurm,
                           template = "batchtools.slurm.tmpl",
                           resources = list(walltime = "10:00:00", ntasks = 1, 
                                            ncpus = 2, memory = 1000)
                       )
                   )
               )),
    tar_rds(save_db, saveRDS(t_set, file = file.path("output", "sets", glue("madin_{t}_{agg}.rds", t = trait, 
                                                                            agg = ifelse(g_agg, "genus", "species")))))
)


list(
    get_path,
    load_db, 
    create_db
)
