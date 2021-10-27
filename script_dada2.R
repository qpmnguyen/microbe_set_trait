library(targets)
library(tarchetypes)
source("R/dada2_func.R")

path_to_data <- function(){
  output <- "~/research/data/ena_gevers_et_al/input_data/"
  return(output)
}


list(
  tar_target(r_data, path_to_data(), format = "file"),
  tar_target(f_trim, {filter_and_trim(r_data)}), 
  tar_target(l_errors, learn_errors(f_trim)),
  tar_target(r_dada, run_dada2(l_errors)),
  tar_target(a_taxonomy, assign_taxonomy(r_dada)),
  tar_rds(save_results, saveRDS(a_taxonomy, file = "data/gevers_dada2.rds"))
)