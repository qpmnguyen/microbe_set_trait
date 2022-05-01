# Performing coverage analysis for each class 
# Quang Nguyen 
# Last updated 2022-05-01  

source(".Rprofile")

# squencing is 16s or wgs an t_category is 
sequencing <- snakemake@params[["sequencing"]]
t_category <- snakemake@params[["t_category"]]


library(TaxaSetsUtils);
library(BiocSet);
library(curatedMetagenomicData);
library(HMP16SData);
library(tidyverse);
library(phyloseq);
library(here);
here::i_am("notebooks/coverage.ipynb")
source(here("R", "generate_scores_func.R"))
source(here("R", "functions_coverage.R"))
options(getClass.msg = FALSE) 


if (sequencing == "16s"){
    data <- HMP16SData::V35() %>% as_phyloseq()
}


