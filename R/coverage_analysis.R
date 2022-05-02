# Performing coverage analysis for each class 
# Quang Nguyen 
# Last updated 2022-05-01  

source(".Rprofile")

# squencing is 16s or wgs and category is one of the trait categories 
sequencing <- snakemake@params[["sequencing"]]
category <- snakemake@params[["category"]]
set_path <- snakemake@input[[1]]

library(TaxaSetsUtils)
library(BiocSet)
library(curatedMetagenomicData)
library(HMP16SData)
library(tidyverse)
library(phyloseq)
library(here)
here::i_am("notebooks/coverage.ipynb")
source(here("R", "generate_scores_func.R"))
source(here("R", "functions_coverage.R"))
options(getClass.msg = FALSE) 


if (sequencing == "16s"){
    data <- HMP16SData::V35() %>% as_phyloseq()
    data <- mapid(obj = data, type = "name") 
    sample_data(data) <- sample_data(data) %>% as(., "data.frame") %>% 
        unite(body_site, c(HMP_BODY_SITE, HMP_BODY_SUBSITE), sep = ":")
    body_sites <- sample_data(data) %>% as(., "data.frame") %>% 
        pull(body_site) %>% unique()
} else if (sequencing == "wgs"){
    data <- sampleMetadata %>% filter(study_name == "HMP_2012") %>% 
        returnSamples(dataType = "relative_abundance", rownames = "NCBI") 
    data <- mia::makePhyloseqFromTreeSummarizedExperiment(data, abund_values = "relative_abundance")
    body_sites <- sampleMetadata %>% filter(study_name == "HMP_2012") %>% 
        pull(body_site) %>% unique()
}

# load sets 
q_sets <- readRDS(set_path)


# initialize results 
results <- vector(mode = "list", length = length(body_sites))
for (i in seq_along(body_sites)){
    print(body_sites[i])
    p_sites <- sample_data(data)$body_site
    subset_physeq <- prune_samples(p_sites == body_sites[i], data) %>% 
        filter_taxa(function(x) sum(x > 0) > 1, TRUE)
    if (sequencing == "wgs"){
        subset_set <- q_sets %>% filter_elementset(element %in% taxa_names(subset_physeq))
    } else if (sequencing == "16s"){
        subset_set <- recode_sets(physeq = subset_physeq, q_sets)
    }
    
    results[[i]] <- eval_func(subset_physeq, q_sets = subset_set, sites = 
                                  body_sites[i], t_class = category)
}

output <- do.call(bind_rows, results)

saveRDS(output, file = snakemake@output[[1]])





