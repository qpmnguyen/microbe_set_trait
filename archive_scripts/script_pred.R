library(targets)
library(tarchetypes)
library(glue)
library(purrr)
library(future)
library(future.callr)
library(future.batchtools)
source("R/functions_pred.R")
tidymodels_prefer()
#plan(multisession)

tar_option_set(workspace_on_error = TRUE)

plan(batchtools_slurm, template = "batchtools.slurm.tmpl")

# evaluation settings 
eval_settings <- cross_df(list(
    type = c("16s", "wgs"),
    benchmark = c("pathway", "traits")
))

# function to load the correct data based on benchmark and type 
load_all_df <- function(class, benchmark){
    if (class == "16s"){
        physeq <- get_gevers()
        if (benchmark == "pathway"){
            df <- get_picrust2(metadata = get_metadata(physeq, class = class), type = benchmark)
        } else if (benchmark == "traits"){
            df <- fit_gsva(physeq, t_rank = "genus")
        }
    } else if (class == "wgs"){
        if (benchmark != "traits"){
            df <- get_ihmp(type = benchmark)
        } else {
            physeq <- get_ihmp(type = benchmark)
            df <- fit_gsva(physeq, t_rank = "species")
        }
    }
    return(df)
}

benchmark <- tar_map(unlist = FALSE, values = eval_settings, 
    tar_target(load_df, load_all_df(class = type, benchmark = benchmark)),
    tar_target(init_wkflow, define_modelflow(data = load_df, n_threads = 5)), 
    tar_target(fit, complete_fit(df = load_df, grid_size = 50, init_wkflow = init_wkflow, 
                                 annotate = paste(class, benchmark, sep = "_"), n_threads = 5), 
               resources = tar_resources(
                   future = tar_resources_future(
                       plan = tweak(
                           batchtools_slurm,
                           template = "batchtools.slurm.tmpl",
                           resources = list(walltime = "10:00:00", ntasks = 1, 
                                            ncpus = 5, memory = 2000)
                       )
                   )
    ))
)

benchmark_combine <- tar_combine(combo_results, benchmark[[3]], 
                                command = dplyr::bind_rows(!!!.x))
save_benchmark <- tar_rds(file_bench, saveRDS(combo_results, 
                                              file = file.path("output","pred", "hmp_wgs_16s_pred.rds")))

list(
    benchmark, 
    benchmark_combine,
    save_benchmark
)

