# This script attempts to generate scores for enrichment sets for each data set using the CBEA package 
# Quang Nguyen 
# Last updated 2022-04-28  
library(here)
library(phyloseq)
library(TaxaSetsUtils)
library(CBEA)
library(tidyverse)
library(speedyseq)
library(taxizedb)
here::i_am("R/generate_scores.R")

physeq <- readRDS(file = here("data","pred_relabun_crc_16s_physeq.rds"))
tax_table(physeq) <- tax_table(physeq) %>% as.data.frame() %>% select(-species) %>% as.matrix() %>% 
test <- TaxaSetsUtils::mapid(physeq, type = "name")
df <- tax_table(physeq) %>% as.data.frame()


g_unq <- g %>% unique() %>% str_trim()
results <- map(g_unq, ~{
    if (is.na(.x)){
        return(NA_character_)
    } else {
        name2taxid(.x, out_type = "summary")
    }
})


src <- src_ncbi()
cmd <- glue("SELECT tax_id")

sql_collect(src, "SELECT * FROM names limit 5")

