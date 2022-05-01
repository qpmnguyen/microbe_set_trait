# This script attempts to generate scores for enrichment sets for each data set using the CBEA package 
# Quang Nguyen 
# Last updated 2022-04-28  

renv::activate()

library(here)
library(phyloseq)
library(TaxaSetsUtils)
library(CBEA)
library(tidyverse)
requireNamespace("speedyseq", quietly = TRUE)
library(BiocSet)
library(glue)
library(mia)
here::i_am("R/generate_scores.R")
source(here("R", "generate_scores_func.R"))

# some sample arguments  
# physeq <- readRDS(file = here("data","pred_relabun_crc_wgs_tse.rds"))
# physeq <- readRDS(file = here("data", "pred_relabun_crc_16s_physeq.rds"))
# condition <- 'crc'
# sequencing <- 'wgs'


physeq_path <- snakemake@input[["physeq"]]
condition <- snakemake@params[["condition"]]
sequencing <- snakemake@params[["sequencing"]]


physeq <- readRDS(file = physeq_path)


if (sequencing == "16s"){
    # remove species level identification 
    # trim whitespaces 
    taxtab <- tax_table(physeq) %>% as.data.frame() %>% 
        select(-species) %>% as.matrix()
    taxtab_trim <- apply(taxtab, 2, str_trim) 
    rownames(taxtab_trim) <- rownames(taxtab)
    tax_table(physeq) <- taxtab_trim 
    # assign ncbiids as columns   
    physeq <- TaxaSetsUtils::mapid(physeq, type = "name")
    
    physeq <- speedyseq::transform_sample_counts(physeq, function(x) x + 10e-5)
    physeq <- speedyseq::transform_sample_counts(physeq, function(x) x/sum(x))
    
    joint_sets <- load_joint_sets(tax_level = "genus")
    asv_sets <- recode_sets(physeq = physeq, joint_set = joint_sets)
    
    tse <- mia::makeTreeSummarizedExperimentFromPhyloseq(physeq)
    
    
    scores <- CBEA::cbea(obj = tse, set = asv_sets, parametric = FALSE, 
                         output = "raw", 
                         abund_values = "counts", adj = FALSE)
    # filter by metadata   
    metadata <- read.csv(file = here("data", glue("pred_picrust2_{condition}_metadata.csv", condition = condition)), 
                         row.names = 1)
    export <- tidy(scores)
    export <- export %>% filter(sample_ids %in% rownames(metadata))
    
} else if (sequencing == "wgs"){
    physeq <- transformSamples(physeq, abund_values = "relative_abundance", 
                              method = "relabundance", 
                              pseudocount = 10e-5)
    joint_sets <- load_joint_sets(tax_level = "species")
    joint_sets <- joint_sets %>% filter_elementset(element %in% rownames(physeq))
    scores <- CBEA::cbea(obj = physeq, set = joint_sets, 
                         output = "raw", parametric = FALSE, 
                         abund_values = "relabundance", adj = FALSE)
    metadata <- read.csv(file = here("data", glue("pred_pathway_{condition}_metadata.csv", condition = condition)), 
                         row.names = 1)
    export <- tidy(scores)
    export <- export %>% filter(sample_ids %in% rownames(metadata))
}

write.csv(export, file = snakemake@output[[1]])




