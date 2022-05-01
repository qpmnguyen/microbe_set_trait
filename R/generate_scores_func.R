# Supplementary functions for generate_scores.R file
# Quang Nguyen
# Last updated 2022-04-29
library(BiocSet)
library(glue)
library(here)

#' This function load sets from all subclasses 
#' @return A joint set 
load_joint_sets <- function(tax_level){
    set_files <- list.files(here("output", "sets"), 
                            pattern = glue("_{tax_level}.rds", 
                                           tax_level = tax_level), 
                            full.names = TRUE)
    # generate sets using ASV identifiers by mapping with NCBI identifiers
    set_list <- purrr::map(set_files, readRDS)
    names(set_list) <- set_files %>% str_extract(glue("(?<=\\/set_).*(?=_{tax_level})", tax_level = tax_level))
    
    # for some reason mutate_elementset is not working for BiocSet 1.8.0
    set_list <- imap(set_list, ~{
        list_obj <- as.list(.x)
        names(list_obj) <- paste(names(list_obj), .y, sep=";")
        set_obj <- BiocSet(list_obj)
        return(set_obj)
    })    
    joint_set <- Reduce(generics::union, set_list)
    return(joint_set)
}



#' This function recodes sets with NCBI ids to the ASV identifiers for a given data set
#' @param physeq The phyloseq object 
#' @param set The BiocSet object
#' @return Another BiocSet tailored to the data.frame
recode_sets <- function(physeq, joint_set){
    asv_ncbi <- tax_table(physeq) %>% as.data.frame() %>% 
        rownames_to_column(var = "asv_id") %>% 
        select(asv_id, ncbi) %>% as_tibble() %>% 
        drop_na(ncbi)
    
    ncbi_ids <- asv_ncbi %>% pull(ncbi)
    reduced_set <- joint_set %>% filter_elementset(element %in% ncbi_ids)
    eset <- es_elementset(reduced_set) %>% as_tibble() %>% 
        dplyr::rename("ncbi" = "element")
    asv_set <- inner_join(asv_ncbi, eset)
    asv_list <- map(unique(asv_set$set), ~{
        asv_set %>% filter(set == .x) %>% pull(asv_id)
    })
    names(asv_list) <- unique(asv_set$set)
    out_set <- BiocSet(asv_list)
    return(out_set)
}