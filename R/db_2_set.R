# Convert databases into sets 
# Last Updated 2022-04-21
# Quang Nguyen  

library(here)
library(tidyverse)
library(BiocSet)
library(rlang)
library(taxizedb)
here::i_am("db_2_set.R")

df <- snakemake@input[["database"]]
class <- snakemake@input[["class"]]
level <- snakemake@input[["level"]]

df <- readRDS(file = here("output", "databases", "db_merged.rds"))
class <- "metabolism"
level <- "species"

if (level == "genus"){
    df <- df %>% mutate(genus_tax_id = map_chr(species_tax_id, ~{
        class <- classification(.x)
        if (is.na(class[[1]])){
            return(NA_character_)
        } else {
            gid <- class[[1]] %>% filter(rank == "genus") %>% pull(id)
            if (length(gid) == 0){
                return(NA_character_)
            } else {
                return(gid)
            }
        }
    }))
    
}




if (class %in% c("metabolism", "gram_stain", "sporulation", "motility", "cell_shape")){
    df <- df %>% select(species_tax_id, {{ class }}) %>% 
        drop_na({{ class }}) 
    if (level == "species"){
        nested_df <- group_by(!!sym(class)) %>% nest(ids = species_tax_id) %>% 
            drop_na({{ class }}) 
        list_ids <- nested_df %>% pull(ids)
        list_ids <- map(list_ids, ~{.x %>% pull(species_tax_id)})
        names(list_ids) <- nested_df %>% pull(class)
        set <- BiocSet::BiocSet(list_ids)
    } else if (level == "genus"){
        classif <- taxizedb::classification(df$species_tax_id)
        
        
        
        new_df <- tibble(
            species_tax_id = names(classif),
            genus_tax_id = genus_id
        )
        g_df <- left_join(df, new_df, by = "species_tax_id")
        g_nested <- g_df %>% group_by(genus_tax_id, !!sym(class)) %>% 
            count()
    }

}



unq <- as.vector(na.omit(unique(df %>% pull(class))))

df %>% group_by(!!class) %>% nest(species_tax_id)


