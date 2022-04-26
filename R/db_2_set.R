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
class <- "metabolism"
level <- "genus"

# Let's extr
if (level == "genus"){
    src <- taxizedb::src_ncbi()
    ncbi_nspec <- sql_collect(src, "SELECT COUNT(*) FROM nodes WHERE rank = 'species'") %>%
        pull()
    db_nspec <- nrow(df)
}


if (class %in% c("metabolism", "gram_stain", "sporulation", 
                 "motility", "cell_shape")){
    df_reduced <- df %>% select(species_tax_id, genus_tax_id, {{ class }}) %>% 
        drop_na({{ class }}) 
    if (level == "species"){
        nested_df <- df_reduced %>% group_by(!!sym(class)) %>% 
            nest(ids = species_tax_id) %>% 
            drop_na({{ class }}) 
        list_ids <- nested_df %>% pull(ids)
        list_ids <- map(list_ids, ~{.x %>% pull(species_tax_id)})
        names(list_ids) <- nested_df %>% pull(class)
        set <- BiocSet::BiocSet(list_ids)
    } else if (level == "genus"){
        
        df_nest <- df_reduced %>% group_by(genus_tax_id, !!sym(class)) %>% count() %>%
            group_by(genus_tax_id) %>% nest(traits = c(metabolism, n))
        start <- Sys.time()
        df_nest <- df_nest %>% mutate(p_vals = test_genus(genus_id = genus_tax_id, full_db = df, 
                                                          db_nspec = db_nspec, ncbi_nspec = ncbi_nspec))
        end <- Sys.time()
        print(paste("Total time is", end - start))
        # adjust for p-values and then filter at 0.05 threshold
        df_nest <- df_nest %>% mutate(p_vals = p.adjust(p_vals, method = "BH")) %>% 
            filter(!is.na(p_vals) | p_vals > 0.05)
        # get species per genera across the database 
        df_genus <- df %>% group_by(genus_tax_id) %>% drop_na(genus_tax_id) %>% 
            count() 
        
        joint_df <- inner_join(df_genus, df_nest, by = "genus_tax_id")
    }

}
saveRDS(set, file = snakemake@output[[1]])



