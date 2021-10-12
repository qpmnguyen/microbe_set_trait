# scripts to pre-process and match the metacyc and 
# taxonomic databases 
library(tidyverse)
library(phyloseq)
library(BiocSet)
library(taxizedb)
library(reticulate)
library(glue)

if (length(tdb_cache$list()) == 0){
  taxizedb::db_download_ncbi(verbose = FALSE, overwrite = FALSE)
}

# call ncbi ids 
# TODO: Modfiy to just return the entire vector 
call_id <- function(ids){
  ncbi <- map(ids, ~{
    res <- name2taxid(.x, out_type = "summary")
    ids <- map_chr(res$id, ~{
      classif <- classification(.x)
      superkingdom <- classif[[1]] |> filter(rank == "superkingdom") |>
        pull(name)
      if (superkingdom == "Bacteria"){
        return(.x)
      } else {
        return(NA_character_)
      }
    })
    ids <- as.vector(na.omit(ids))
    if (length(ids) >= 2){
      message("Ambiguous identifiers")
    }
    return(ids)
  })
  return(ncbi)
}


#' @param trait_vec Vector of unique traits 
#' @param db Database where species/taxa as rows and traits as columns.
#' @param filt_term Term to filter database on 
#' @param id_col Name of id column indicating NCBI identifiers 
#' @param genus_agg What are the direct levels 
#' @param threshold Threshold for when a trait is considered
retr_sets <- function(db, trait_vec, filt_term, 
                      id_col, genus_agg = FALSE, threshold = 0.95){
  # special processing if filtered term is carbon_substrates, pathways
  set_list <- map(trait_vec, ~{
    if (filt_term %in% c("carbon_substrates", "pathways")){
      db_inter <- db |> 
        rowwise() |> 
        mutate(t_bool = str_detect(string = !!sym(filt_term), 
                                   pattern = .x)) |> 
        ungroup() 
    } else{
      db_inter <- db |> 
        mutate(t_bool = ifelse(!!sym(filt_term) == .x, TRUE, FALSE))
    }
    if (genus_agg == TRUE){
      db_inter |> group_by(genus) |> 
        summarise(t_count = sum(t_bool, na.rm = TRUE), n = n()) |> 
        mutate(prop = t_count/n) |> 
        filter(prop >= threshold) |> 
        mutate(ncbi_id = call_id(genus)) |> 
        pull(ncbi_id) |>
        flatten_chr()
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
process_db <- function(){
  trait_db <- read.table(file = "data/condensed_species_NCBI.txt", sep = ",", header = TRUE)
  metadata <- colnames(trait_db)[1:8]
  traits <- c("metabolism", "pathways", "carbon_substrates", "sporulation")
  # restrict databases to certain types of traits 
  trait_db <- trait_db |> as_tibble() |>
    select(all_of(metadata), all_of(traits)) |> 
    dplyr::filter(superkingdom == "Bacteria") 
  
  return(trait_db)
}

get_trait_list <- function(trait_db){
  traits <- c("metabolism", "pathways", "carbon_substrates", "sporulation")
  trait_list <- map(traits, ~{
    t_vec <- unique(trait_db |> pull(!!sym(.x)))
    if (.x %in% c("pathways", "carbon_substrates")){
      t_list <- map(t_vec, ~ str_split(.x, ",")) |> flatten() |> flatten_chr()
      t_vec <- map_chr(t_list, ~ str_trim(.x)) |> unique() |> na.omit() |> as.vector()
    } else {
      t_vec <- na.omit(t_vec)
    }
    return(t_vec)
  })
  return(trait_list)
}


create_sets <- function(trait_list, trait_db, agg = FALSE, threshold = 0.95){
  traits <- c("metabolism", "pathways", "carbon_substrates", "sporulation")
  trait_sets <- imap(trait_list, ~{
    retr_sets(db = trait_db, trait_vec = .x, filt_term = traits[.y], 
              id_col = "species_tax_id", genus_agg = agg, threshold = thresh)
  })
}

# do not run 
process_metacyc <- function(){
  # please make sure pathway tools is running in the background
  # The following commands did not work because it forces
  # the current R console to be in the pathway-tools program
  # with no easy way to exit
  # pthway_tools_exec <- "~/pathway-tools"
  # system(glue("{dir}/pathway-tools -lisp -python-local-only", 
  #            dir = pthway_tools_exec))
  Sys.setenv(RETICULATE_PYTHON = "~/miniconda3/envs/microbe_set_trait/bin/python3")
  reticulate::source_python("python/metacyc.py")
  # get_hierarchy and get_instances are python functions
  subclass <- get_hierarchy()  
  inst <- map(subclass, ~{
    get_instances(.x)
  })
  names(inst) <- subclass
  return(inst)
}



