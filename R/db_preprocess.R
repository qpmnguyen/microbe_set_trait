# Functions to process the finalized database into BiocSets for both genus and species level annotation
# Last modified 2022-04-26
# Quang Nguyen 

library(tidyverse)
library(BiocSet)
library(taxizedb)
library(glue)

#' Function to assess whether a genus is underrepresented in this data frame
genus_assess <- function(df_reduced, full_db, category, db_nspec, ncbi_nspec, filt_thresh){
    # obtain nesting data frame and test using the hypergeometric test for under representation
    # nest traits within each genus. If trait is a vector list, unnest it first
    df_nest <- df_reduced %>% unnest(trait) %>%
        group_by(genus_tax_id) %>% nest(trait_set = c(species_tax_id, trait)) %>% ungroup() 
    
    # obtain p-values by performing the hypergeometric test 
    df_nest <- df_nest %>% mutate(p_vals = test_genus(genus_id = genus_tax_id, df_reduced = df_reduced, 
                                                      db_nspec = db_nspec, ncbi_nspec = ncbi_nspec))
    # adjust for p-values and filter which ones are not significant at 0.05 level 
    df_nest <- df_nest %>% mutate(p_vals = p.adjust(p_vals, method = "BH")) %>% 
        filter(!is.na(p_vals)) %>% filter(p_vals > 0.05)
    
    # df_genus represent the total number of species with an annotated trait of the same 
    # category per genus
    df_genus <- full_db %>% drop_na(genus_tax_id) %>%
        drop_na(!!sym(category)) %>% group_by(genus_tax_id) %>%
        count() 
    
    # join them in to compute proportions
    joint_df <- inner_join(df_genus, df_nest, by = "genus_tax_id") %>% ungroup()
    
    # count number of species having the trait for each trait 
    final_df <- joint_df %>% mutate(trait_set = map(trait_set, ~{
        .x %>% group_by(trait) %>% count() %>% 
            rename("nspec" = "n") %>%
            ungroup()
    })) %>% unnest(trait_set)
    
    final_df <- final_df %>%
        mutate(prop = nspec/n) %>% 
        filter(prop >= filt_thresh) %>% 
        select(-c(n, nspec, p_vals, prop)) %>% 
        group_by(trait) %>% summarize(ids = list(genus_tax_id))
    # for each trait that has the proportion at least 95% of total species having
    # at least one trait within the category 
    list_ids <- final_df %>% pull(ids)
    names(list_ids) <- final_df %>% pull(trait)
    str_replace_all(names(list_ids), pattern = " ", replacement = "_")
    return(list_ids)
}

#' n_spec is not part of the function since it is an expensive operation that
#' should be done only once
#' @param genus_id NCBI identifier for genus as string or numeric values. If
#'     string, convert to numeric
#' @param db_nspec Total number of species in the trait database (nrow of the database), 
#'     considered to be equivalent to the paramter k in a hypergeometric distribution
#' @param db_ngenus Total number of species assigned to the genus within the database, 
#'     considered to be equivalent to the parameter q in a hypergeometric distribution. 
#' @param ncbi_nspec Total number of species in NCBI. Considered to be m + n parameter (or nn)
#'     in a hypergeometric distribution. 
test_genus <- function(genus_id, df_reduced, db_nspec, ncbi_nspec) {
    # check vector 
    if (!is.vector(genus_id)) {
        genus_id <- as.vector(genus_id)
    }
    genus_id <- map_dbl(genus_id, as.numeric)
    # obtain p-values
    p_vals <- map_dbl(genus_id, function(x){
        if (is.na(x)){
            return(NA_real_)
        } else {
            # among the species with annotated trait, how many species are assigned
            # to the genus
            db_ngenus <- df_reduced %>% filter(genus_tax_id == x) %>% nrow()
            cmd <- paste("SELECT tax_id, level, ancestor FROM hierarchy WHERE ancestor =", x)
            tbl <- sql_collect(src, cmd)
            ncbi_ngenus <- tbl %>% mutate(rank = taxid2rank(tax_id)) %>% filter(rank == "species") %>%
                nrow()
            #print(ncbi_ngenus)
            #print(db_ngenus)
            cum_prob <- phyper(
                q = db_ngenus, 
                m = ncbi_ngenus, 
                n = ncbi_nspec - ncbi_ngenus,
                k = db_nspec, 
                lower.tail = TRUE
            )
            return(cum_prob)
        }
    })
    names(p_vals) <- genus_id
    return(p_vals)
}

# this function process a list of pvalues and returns NAs if
# things go wrong
proc_gtest <- function(pvals) {
    if (length(pvals) == 0) {
        return(NA_character_)
    }
    sig <- which(pvals < 0.05)
    if (length(sig) == 0) {
        return(NA_character_)
    } else {
        return(names(pvals)[sig])
    }
}



# PROCESS NCBIID RESULTS ####  
# do not run if trying to perform this reproducibly
# this is due to the fact that it requires an external software
# Please visit pathwaytools page for more information.
process_metacyc <- function() {
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
    inst <- map(subclass, ~ {
        get_instances(.x)
    })
    names(inst) <- subclass
    return(inst)
}
