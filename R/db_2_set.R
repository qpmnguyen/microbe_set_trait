# Convert databases into sets 
# Last Updated 2022-04-25
# Quang Nguyen  

library(here)
library(tidyverse)
library(BiocSet)
library(rlang)
library(taxizedb)
library(glue)
here::i_am("R/db_2_set.R")
source(here("R", "db_preprocess.R"))

df <- snakemake@input[["database"]]
class <- snakemake@input[["class"]]
level <- snakemake@input[["level"]]

df <- readRDS(file = here("output", "databases", "db_merged.rds"))
class <- "substrate"
level <- "genus"

# Obtain some constants 
if (level == "genus"){
    src <- taxizedb::src_ncbi()
    ncbi_nspec <- sql_collect(src, "SELECT COUNT(*) FROM nodes WHERE rank = 'species'") %>%
        pull()
    db_nspec <- nrow(df)
}


if (class %in% c("metabolism", "gram_stain", "sporulation", 
                 "motility", "cell_shape")){
    df_reduced <- df %>% select(species_tax_id, genus_tax_id, {{ class }}) %>% 
        rename("trait" = class) %>%
        drop_na(trait) 
    if (level == "species"){
        # for species we simply extract the species associated with each trait 
        nested_df <- df_reduced %>% group_by(trait) %>% 
            nest(ids = species_tax_id) %>% 
            drop_na(trait) 
        list_ids <- nested_df %>% pull(ids)
        list_ids <- map(list_ids, ~{.x %>% pull(species_tax_id)})
        names(list_ids) <- nested_df %>% pull(trait)
        set <- BiocSet::BiocSet(list_ids)
    } else if (level == "genus"){
        # here, first, we nest traits within genera and test for under-representation
        df_nest <- df_reduced %>% group_by(genus_tax_id, trait) %>% count() %>%
            group_by(genus_tax_id) %>% nest(trait_set = c(trait, n))
        start <- Sys.time()
        df_nest <- df_nest %>% mutate(p_vals = test_genus(genus_id = genus_tax_id, full_db = df, 
                                                          db_nspec = db_nspec, ncbi_nspec = ncbi_nspec))
        end <- Sys.time()
        print(paste("Total time is", format(end - start)))
        # adjust for p-values and then filter at 0.05 threshold
        df_nest <- df_nest %>% mutate(p_vals = p.adjust(p_vals, method = "BH")) %>% 
            filter(!is.na(p_vals) | p_vals > 0.05)
        # get species per genera across the database 
        df_genus <- df %>% group_by(genus_tax_id) %>% drop_na(genus_tax_id) %>% 
            count() 
        # then we test whether at least 95% of species of that genera in the ENTIRE DATABASE has the trait
        # this means we definitely count NAs to calculate proportions. 
        joint_df <- inner_join(df_genus, df_nest, by = "genus_tax_id")
        
        joint_df <- joint_df %>% mutate(traits = map2_chr(traits, n, ~{
            true_prop <- .x %>% mutate(prop = n / .y) %>% 
                filter(prop >= 0.95)
            if (nrow(true_prop) >= 2 | nrow(true_prop) == 0){
                return(NA_character_)
            } else {
                return(true_prop %>% pull(class))
            }
        })) 
        
        joint_df <- joint_df %>% select(-c(n, p_vals)) %>% 
            na.omit() %>% group_by(traits) %>% nest(ids = genus_tax_id) 
        list_ids <- joint_df %>% pull(ids) %>% map(., ~{pull(.x, genus_tax_id)})
        names(list_ids) <- joint_df %>% pull(traits)
        set <- BiocSet(list_ids)
    }
} else if (class %in% c("pathways", "substrate")){
    df_reduced <- df %>% select(species_tax_id, genus_tax_id, !!sym(class)) %>%
        rename("trait" = class) %>%
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
        set <- BiocSet(list_ids)
    } else if (level == "genus"){
        joint_df <- genus_assess(df_reduced = df_reduced, full_db = df)
        
    }
}
saveRDS(set, file = snakemake@output[[1]])












