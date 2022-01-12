library(targets)
library(tarchetypes)
library(glue)
library(future)
library(future.callr)

source("R/utils.R")
source("R/functions_coverage.R")

plan(callr)
eval_settings <- function(analysis){
    if (analysis == "hmp_16s"){
        physeq <- retrieve_hmp_16s()
    } else if (analysis == "hmp_wgs") {
        physeq <- retrieve_hmp_wgs()
    }
    settings <- cross_df(list(
        site = unique(sample_data(physeq) %>% pull(body_site)),
        type = c("carbon_substrates", "cell_shape", "gram_stain", "motility", 
                 "pathways", "range_salinity", "range_tmp", "sporulation")
    ))
    return(settings)
}

physeq_hmp_16s <- tar_target(hmp_16s, retrieve_hmp_16s())
mapped_hmp_16s <- tar_map(unlist = FALSE, values = eval_settings("hmp_16s"), 
        tar_target(load_hmp_16s, {
            physeq <- subset_samples(hmp_16s, body_site == site)
            # filter so that each taxa is at at least one sample 
            physeq <- filter_taxa(physeq, function(x) sum(x > 0) >= 1)
            physeq <- attach_ncbi_std(physeq, t_rank = "GENUS")
            return(physeq)
        }),
        tar_target(load_sets, {
            readRDS(file = file.path("output", "sets", glue("madin_{type}_genus.rds", type = type)))
        }),
        tar_target(hmp_16s_cov, {
            calculate_coverage(physeq = load_hmp_16s, sets = load_sets, type = type, site = site)
        })
)
combined_hmp_16s <- tar_combine(comb_hmp_16s, mapped_hmp_16s[[3]], 
                                command = dplyr::bind_rows(!!!.x))
save_hmp_16s <- tar_rds(file_hmp_16s, saveRDS(comb_hmp_16s, 
                                              file = file.path("output","coverage", "hmp_16s_cov.rds")))

physeq_hmp_wgs <- tar_target(hmp_wgs, retrieve_hmp_wgs())
mapped_hmp_wgs <- tar_map(unlist = FALSE, values = eval_settings("hmp_wgs"),
        tar_target(load_hmp_wgs, {
            physeq <- subset_samples(hmp_wgs, body_site == site)
            # filter so that each taxa is at at least one sample 
            physeq <- filter_taxa(physeq, function(x) sum(x > 0) >= 1)
            physeq <- attach_ncbi_metaphlan(physeq)
            return(physeq)
        }),
        tar_target(load_sets_2, {
            readRDS(file = file.path("output", "sets", glue("madin_{type}_species.rds", type = type)))
        }),
        tar_target(hmp_wgs_cov, {
            calculate_coverage(physeq = load_hmp_wgs, sets = load_sets_2, type = type, site = site)
        })
)
combined_hmp_wgs <- tar_combine(comb_hmp_wgs, mapped_hmp_wgs[[3]], 
                                command = dplyr::bind_rows(!!!.x))
save_hmp_wgs <- tar_rds(file_hmp_wgs, saveRDS(comb_hmp_wgs, 
                                              file = file.path("output","coverage", "hmp_16s_cov.rds")))



# list(

# )
list(
    physeq_hmp_wgs,
    mapped_hmp_wgs,
    combined_hmp_wgs,
    save_hmp_wgs,
    physeq_hmp_16s,
    mapped_hmp_16s,
    combined_hmp_16s,
    save_hmp_16s
)




