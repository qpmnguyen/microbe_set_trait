library(HMP16SData)
library(curatedMetagenomicData)
library(BiocSet)
library(magrittr)
library(tibble)
library(dplyr)
library(tidyr)
library(purrr)
source("R/db_preprocess.R")
source("R/utils.R")

requireNamespace("speedyseq", quietly = TRUE)
# DATA RETRIEVAL ####
retrieve_hmp_16s <- function(){
    df <- list(V13 = V13() , V35 = V35())
    phy_list <- map(df, ~{
        physeq <- .x %>% subset(select = VISITNO == 1) %>% as_phyloseq()
        physeq <- prune_samples(sample_sums(physeq) >= 1000, physeq)
        #physeq <- filter_taxa(physeq, function(x) sum(x > 0) > 0.05 * length(x), TRUE)
        phyloseq(
            otu_table = otu_table(physeq),
            sample_data = sample_data(physeq),
            tax_table = tax_table(physeq)
        )
    })
    physeq <- merge_phyloseq(phy_list[[1]], phy_list[[2]])
    sample_data(physeq) <- sample_data(physeq) %>% as(.,"data.frame") %>% 
        unite(body_site, c(HMP_BODY_SITE, HMP_BODY_SUBSITE), sep = ":")
    rm(phy_list)
    rm(df)
    gc()
    return(physeq)
}

retrieve_hmp_wgs <- function(){
    physeq <- suppressMessages(curatedMetagenomicData(pattern = "HMP_2012.relative_abundance", dryrun = FALSE)) %>% 
        mergeData()
    physeq <- suppressMessages(mia::makePhyloseqFromTreeSummarizedExperiment(physeq, 
                                                                           abund_values = "relative_abundance"))
    physeq <- phyloseq(
        otu_table = otu_table(physeq),
        sample_data = sample_data(physeq),
        tax_table = tax_table(physeq)
    )
    sample_data(physeq) <- sample_data(physeq) %>% as(.,"data.frame") %>% 
        unite(col = "body_site", c(body_site, body_subsite), sep = ":")
    return(physeq)
}

# DATA PROCESSING ####
#' @title Attaching NCBIIDs to a phyloseq data set from a certain rank 
attach_ncbi_std <- function(physeq, t_rank){
    c_ranks <- rank_names(physeq)
    t_rank <- match.arg(t_rank, c_ranks)
    c_ranks <- c_ranks[seq_len(which(c_ranks == t_rank))]
    physeq <- speedyseq::tax_glom(physeq = physeq, taxrank = t_rank)
    taxtab <- tax_table(physeq)
    taxtab <- as(taxtab, "matrix")[,c_ranks]
    
    # this returns a vector 
    ncbiids <- map_chr(seq_len(nrow(taxtab)), ~{
        calls <- call_genus_id(taxtab[.x,])
        if (is.null(calls) | length(calls) == 0){
            return(NA_character_)
        } else {
            return(calls)
        }
    })
    taxtab <- cbind(taxtab, ncbiids)
    tax_table(physeq) <- taxtab
    return(physeq)
}

attach_ncbi_metaphlan <- function(physeq){
    taxtab <- as(tax_table(physeq), "data.frame")
    query_db <- readRDS(file = file.path("mpa_marker.rds"))
    query_db <- query_db %>% as_tibble() %>% select(ncbiID, Species)
    colnames(taxtab) <- stringr::str_to_title(colnames(taxtab))
    taxtab <- taxtab %>% rownames_to_column(var = "full_name") %>% 
        mutate(across(c(Phylum, Class, Order, Family, Genus, Species), 
                      ~str_replace_all(.x, pattern = " ", replacement = "_"))) %>%
        left_join(query_db, by = "Species")  %>%
        as.data.frame() %>% 
        column_to_rownames(var = "full_name") %>% 
        mutate(ncbiID = str_trim(ncbiID)) %>% 
        rename(ncbiids = ncbiID) %>%
        as.matrix()
    tax_table(physeq) <- taxtab
    return(physeq)
}

# CALCULATE COVERAGE ####
#' @title Calculate the desired coverage counts for the body sites  
calculate_coverage <- function(physeq, sets, type, site){
    missing_ncbi <- ntaxa(physeq) - physeq %>% 
        subset_taxa(!is.na(ncbiids)) %>% ntaxa()
    physeq <- physeq %>% subset_taxa(!is.na(ncbiids))
    df_ids <- tax_table(physeq) %>% as(., "data.frame") %>% pull("ncbiids")
    sets_lab <- es_set(sets) %>% pull(set) %>% unique()
    set2df <- map(sets_lab,~{
        set_ids <- es_elementset(sets) %>% filter(set == .x) %>% 
            pull(element)
        length(intersect(df_ids, set_ids))/length(df_ids)
    })
    names(set2df) <- sets_lab
    
    
    out <- tibble(
        missing_ncbi = missing_ncbi,
        tot_match = length(intersect(df_ids, es_element(sets) %>% pull(element)))/length(df_ids), 
        by_set = list(as_tibble(set2df) %>% 
                          pivot_longer(cols = everything(), 
                                       names_to = "set_names", values_to = "prop")),
        type = type,
        site = site
    )
    return(out)
}
