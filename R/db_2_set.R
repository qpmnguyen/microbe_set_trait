# Convert databases into sets primary inputs for snakemake
# Last Updated 2022-04-26
# Quang Nguyen
source(".Rprofile")

library(here)
library(tidyverse)
library(BiocSet)
library(rlang)
library(taxizedb)
library(glue)
here::i_am("R/db_2_set.R")
source(here("R", "db_preprocess.R"))

# These are test parameters
#df <- readRDS(file = here("output", "databases", "db_merged.rds"))
#category <- "metabolism"
#level <- "genus"
#filt_thresh <- 0.95

df <- snakemake@input[["database"]]
category <- snakemake@params[["category"]]
level <- snakemake@params[["level"]]
filt_thresh <- snakemake@params[["filt_thresh"]]


# Obtain some constants 
if (level == "genus"){
    src <- taxizedb::src_ncbi()
    ncbi_nspec <- sql_collect(src, "SELECT COUNT(*) FROM nodes WHERE rank = 'species'") %>%
        pull()
}

df <- readRDS(df)

if (category %in% c("metabolism", "gram_stain", "sporulation", 
                 "motility", "cell_shape")){
    # obtain a reduced data set by removing all species with NA traits
    # and then rename into one column (trait)
    df_reduced <- df %>% select(species_tax_id, genus_tax_id, {{ category }}) %>% 
        rename("trait" = category) %>%
        drop_na(trait) 
    
    if (level == "species"){
        # for species we simply extract the species associated with each trait 
        nested_df <- df_reduced %>% group_by(trait) %>% 
            summarise(ids = list(species_tax_id)) %>% 
            drop_na(trait) 
        list_ids <- nested_df %>% pull(ids)
        names(list_ids) <- nested_df %>% pull(trait)
        str_replace_all(names(list_ids), pattern = " ", replacement = "_")
        set <- BiocSet::BiocSet(list_ids)
    } else if (level == "genus"){
        print(filt_thresh)
        db_nspec <- df %>% drop_na(!!sym(category)) %>% drop_na(genus_tax_id) %>% nrow()
        list_ids <- genus_assess(df_reduced = df_reduced, full_db = df, category = category, 
                                 db_nspec = db_nspec, ncbi_nspec = ncbi_nspec, 
                                 filt_thresh = filt_thresh)
        set <- BiocSet(list_ids)
    }
} else if (category %in% c("pathways", "substrate")){
    df_reduced <- df %>% select(species_tax_id, genus_tax_id, !!sym(category)) %>%
        rename("trait" = category) %>%
        mutate(trait = na_if(trait, "NA")) %>%
        drop_na(trait) %>% 
        mutate(trait = map(trait, ~str_split(.x, ",")[[1]]))
    
    trait_list <- df_reduced %>% pull(trait) %>% flatten_chr() %>% unique()
    
    if (level == "species") {
        list_ids <- vector(mode = "list", length = length(trait_list))
        for (i in seq_along(trait_list)){
            list_ids[[i]] <- df_reduced %>% 
                mutate(check = map_lgl(trait, ~{trait_list[i] %in% .x})) %>% 
                filter(check == TRUE) %>% pull(species_tax_id)
        }
        names(list_ids) <- trait_list
        str_replace_all(names(list_ids), pattern = " ", replacement = "_")
        set <- BiocSet(list_ids)
    } else if (level == "genus"){
        db_nspec <- df_reduced %>% drop_na(genus_tax_id) %>% nrow()
        list_ids <- genus_assess(df_reduced = df_reduced, full_db = df, category = category, db_nspec = db_nspec, 
                                 ncbi_nspec = ncbi_nspec, filt_thresh = filt_thresh)
        set <- BiocSet(list_ids)
    }
}
saveRDS(set, file = snakemake@output[[1]])