#' Master script to run all individual pipelines  
library(targets)
library(optparse)

# run the dada2 pipeline to process gevers et al
# change input directory if data is posted elsewhere 
tar_make(script = "script_dada2.R", store = "store_dada2")
tar_make(script = "script_db_preprocess.R", store = "store_db_preprocess")
tar_make(script = "script_prep_data.R", store = "store_prep_data")


tar_visnetwork(script = "script_dada2.R", store = "store_dada2")
