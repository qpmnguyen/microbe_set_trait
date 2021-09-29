library(tidyverse)
library(phyloseq)
library(BiocSet)
library(taxizedb)

if (!file.exists(tdb_cache$list())){
  taxizedb::db_download_ncbi(verbose = FALSE, overwrite = FALSE)
}

#' @param trait_vec Vector of unique traits 
#' @param db Database where species/taxa as rows and traits as columns.
#' @param filt_term Term to filter database on 
#' @param id_col Name of id column indicating NCBI identifiers 
#' @param genus_agg What are the direct levels 
#' @param threshold Threshold for when a trait is considered
retr_sets <- function(db, trait_vec, filt_term, id_col, genus_agg, threshold = 0.95){
  if (filt_term %in% c("carbon_substrates", "pathways")){
    set_list <- map(trait_vec, ~{
      db |> filter(.x %in% !!sym(filt_term)) |> 
        pull(!!sym(id_col)) |> as.character()
    })
  } else {
    set_list <- map(trait_vec, ~{
      db |> filter(!!sym(filt_term) == .x) |> 
        pull(!!sym(id_col)) |> as.character()
    })
  }
  names(set_list) <- trait_vec
  sets <- BiocSet(set_list)
  sets <- sets |> mutate_set(type = filt_term)
  
  return(sets)
}

#' This function returns a processed trait_set 
process_db <- function(){
  trait_db <- read.table(file = "data/condensed_species_NCBI.txt", sep = ",", header = TRUE)
  metadata <- colnames(trait_db)[1:8]
  traits <- c("metabolism", "pathways", "carbon_substrates", "sporulation")
  
  # restrict databases to certain types of traits 
  trait_db <- trait_db |> as_tibble() |>
    select(all_of(metadata), all_of(traits)) |> 
    dplyr::filter(superkingdom == "Bacteria") 
  
  trait_list <- map(traits, ~{
    t_vec <- unique(trait_db |> pull(!!sym(.x)))
    if (.x %in% c("pathways", "carbon_substrates")){
      t_list <- map(t_vec, ~ str_split(.x, ",")) |> flatten() |> flatten_chr()
      t_vec <- map_chr(t_list, ~ str_trim(.x)) |> unique() |> na.omit()
    } else {
      t_vec <- na.omit(t_vec)
    }
    return(t_vec)
  })
  
  trait_sets <- imap(trait_list, ~{
    retr_sets(db = trait_db, trait_vec = .x, filt_term = traits[.y], 
              id_col = "species_tax_id")
  })
  return(trait_sets)
}




