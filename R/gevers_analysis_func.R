library(tidyverse)
library(phyloseq)

process_db <- function(){
  trait_db <- read.table(file = "data/condensed_species_NCBI.txt", sep = ",", header = TRUE)
  metadata <- colnames(trait_db)[1:8]
  
  trait_db <- trait_db %>% as_tibble() %>% 
    select(all_of(metadata), metabolism, pathways, carbon_substrates) %>% 
    dplyr::filter(superkingdom == "Bacteria")
  
  
  
  pathway_sets <- unique(trait_db$pathways)
  
  pathway_sets <- pathway_sets %>% map(~{str_split(.x, ",")}) %>% flatten() %>% flatten_chr() %>% 
    map_chr(~{str_trim(.x)}) %>% unique()
  
  

  
}