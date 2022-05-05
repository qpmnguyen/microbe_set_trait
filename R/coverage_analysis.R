# Performing coverage analysis for each class 
# Quang Nguyen 
# Last updated 2022-05-01  

source(".Rprofile")

# squencing is 16s or wgs and category is one of the trait categories 
dset <- snakemake@params[["dset"]]
category <- snakemake@params[["category"]]
set_path <- snakemake@input[[1]]

library(TaxaSetsUtils)
library(BiocSet)
library(curatedMetagenomicData)
library(HMP16SData)
library(tidyverse)
library(phyloseq)
library(here)
source(here("R", "generate_scores_func.R"))
source(here("R", "functions_coverage.R"))
options(getClass.msg = FALSE) 


if (dset == "hmp_16s"){
    data <- HMP16SData::V35() %>% as_phyloseq()
    data <- mapid(obj = data, type = "name") 
    sample_data(data) <- sample_data(data) %>% as(., "data.frame") %>% 
        unite(body_site, c(HMP_BODY_SITE, HMP_BODY_SUBSITE), sep = ":")
    body_sites <- sample_data(data) %>% as(., "data.frame") %>% 
        pull(body_site) %>% unique()
} else if (dset == "hmp_wgs"){
    data <- sampleMetadata %>% filter(study_name == "HMP_2012") %>% 
        returnSamples(dataType = "relative_abundance", rownames = "NCBI") 
    data <- mia::makePhyloseqFromTreeSummarizedExperiment(data, abund_values = "relative_abundance")
    body_sites <- sampleMetadata %>% filter(study_name == "HMP_2012") %>% 
        pull(body_site) %>% unique()
} else if (dset == "crc_16s"){
    data <- readRDS(file = here("data", "pred_relabun_crc_16s.rds"))
    taxtab <- tax_table(data) %>% as.data.frame() %>% 
        select(-species) %>% as.matrix()
    taxtab_trim <- apply(taxtab, 2, str_trim) 
    rownames(taxtab_trim) <- rownames(taxtab)
    tax_table(data) <- taxtab_trim
    data <- mapid(obj = data, type = "name")
} else if (dset == "crc_wgs") {
    data <- readRDS(file = here("data", "pred_relabun_crc_wgs.rds"))
    data <- mia::makePhyloseqFromTreeSummarizedExperiment(data, abund_values = "relative_abundance")
} else if (dset == "ibd_16s"){
    data <- readRDS(file = here("data", "pred_relabun_ibd_16s.rds"))
    taxtab <- tax_table(data) %>% as.data.frame() %>% 
        select(-species) %>% as.matrix()
    taxtab_trim <- apply(taxtab, 2, str_trim) 
    rownames(taxtab_trim) <- rownames(taxtab)
    tax_table(data) <- taxtab_trim
    data <- mapid(obj = data, type = "name")
} else if (dset == "ibd_wgs"){
    data <- readRDS(file = here("data", "pred_relabun_ibd_wgs.rds"))
    data <- mia::makePhyloseqFromTreeSummarizedExperiment(data, abund_values = "relative_abundance")
}

# load sets 
q_sets <- readRDS(set_path)

if (str_detect(dset, "hmp")){
    # loop through body sites if hmp is being evaluated 
    # initialize results 
    results <- vector(mode = "list", length = length(body_sites))
    for (i in seq_along(body_sites)){
        print(body_sites[i])
        p_sites <- sample_data(data)$body_site
        subset_physeq <- prune_samples(p_sites == body_sites[i], data) %>% 
            filter_taxa(function(x) sum(x > 0) > 1, TRUE)
        subset_physeq <- prune_samples(sample_sums(subset_physeq) != 0, subset_physeq)
        if (str_detect(dset, "wgs")){
            subset_set <- q_sets %>% filter_elementset(element %in% taxa_names(subset_physeq))
        } else if (str_detect(dset, "16s")){
            subset_set <- recode_sets(physeq = subset_physeq, q_sets)
        }
        
        results[[i]] <- eval_func(subset_physeq, q_sets = subset_set, sites = 
                                      body_sites[i], t_class = category)
    }
    output <- do.call(bind_rows, results)
} else {
    subset_physeq <- data %>% filter_taxa(function(x) sum(x > 0) > 1, TRUE)
    subset_physeq <- prune_samples(sample_sums(subset_physeq) != 0, subset_physeq)
    if (str_detect(dset, "wgs")){
        subset_set <- q_sets %>% filter_elementset(element %in% taxa_names(subset_physeq))
    } else if (str_detect(dset, "16s")){
        subset_set <- recode_sets(physeq = subset_physeq, q_sets)
    }
    output <- eval_func(subset_physeq, q_sets = subset_set, sites = 
                  dset, t_class = category)
}

saveRDS(output, file = snakemake@output[[1]])





