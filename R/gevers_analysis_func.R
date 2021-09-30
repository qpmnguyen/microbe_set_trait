library(tidyverse)
library(phyloseq)
library(BiocSet)
library(taxizedb)

if (!file.exists(tdb_cache$list())){
  taxizedb::db_download_ncbi(verbose = FALSE, overwrite = FALSE)
}

# call ncbi ids 
# TODO: Modfiy to just return the entire vector 
call_id <- function(ids){
  ncbi <- map_chr(ids, ~{
    res <- name2taxid(.x, out_type = "summary")
    return(as.character(res[1,2]))
  })
  return(ncbi)
}


#' @param trait_vec Vector of unique traits 
#' @param db Database where species/taxa as rows and traits as columns.
#' @param filt_term Term to filter database on 
#' @param id_col Name of id column indicating NCBI identifiers 
#' @param genus_agg What are the direct levels 
#' @param threshold Threshold for when a trait is considered
retr_sets <- function(db, trait_vec, filt_term, id_col, genus_agg = FALSE, threshold = 0.95){
  # special processing if filtered term is carbon_substrates, pathways
  set_list <- map(trait_vec, ~{
    if (filt_term %in% c("carbon_substrates", "pathways")){
      db_inter <- db |> 
        rowwise() |> 
        mutate(t_bool = str_detect(string = !!sym(filt_term), pattern = .x)) |> 
        ungroup() 
    } else{
      db_inter <- db |> mutate(t_bool = !!sym(filt_term) == .x)
    }
    if (genus_agg == TRUE){
      db_inter |> group_by(genus) |> 
        summarise(t_count = sum(t_bool, na.rm = TRUE), n = n()) |> 
        mutate(prop = t_count/n) |> 
        filter(prop >= threshold) |> 
        mutate(ncbi_id = call_id(genus)) |> 
        pull(ncbi_id) |> as.character()
    } else {
      db_inter |> filter(t_bool == TRUE) |> pull(!!sym(id_col)) |> as.character()
    }
  })
  names(set_list) <- trait_vec
  sets <- BiocSet(set_list)
  sets <- sets |> mutate_set(type = filt_term)
  
  return(sets)
}

#' This function returns a processed trait_set 
process_db <- function(agg = FALSE, thresh = 0.95){
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
              id_col = "species_tax_id", genus_agg = agg, threshold = thresh)
  })
  return(trait_sets)
}

process_metacyc <- function(){
  # change db_path
  db_path <- file.path("databases", "meta_25.1", "data")
  if (length(list.files(db_path)) == 0){
    stop("Need to download metaCyc pathways to db_path")
  }
  
  class_file <- file.path(db_path, "classes.dat")
  classes <- readr::read_file(class_file)
  
  
}

