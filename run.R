#' Master script to run all individual pipelines  
library(targets)

# run the dada2 pipeline to process gevers et al
# change input directory if data is posted elsewhere 
tar_make(script = "script-dada2.R", store = "store_dada2")

