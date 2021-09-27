library(tidyverse)
library(phyloseq)

process_db <- function(){
  trait_db <- read.table(file = "data/condensed_species_NCBI.txt", sep = ",", header = TRUE)
  metadata <- colnames(trait_db)[1:8]
  
  trait_db <- trait_db %>% as_tibble() %>% 
    select(all_of(metadata), metabolism, pathways, carbon_substrates) %>% 
    dplyr::filter(superkingdom == "Bacteria")
  
  
}