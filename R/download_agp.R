# Created March 7th 2022
# Quang Nguyen 
# This script downloads and merges temporary

library(tidyverse)
library(piggyback)
library(here)
here::i_am("R/download_agp.R")


f1 <- tempdir()
f2 <- tempfile()

# this is the metadata file  
pb_download(file = "agp_metadata_20220304.txt", dest = f1, tag = "0.1")

# file manifest from ebi-ena
download.file(url = "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJEB11419&result=read_run&fields=study_accession,sample_accession,run_accession,library_layout,library_strategy,library_source,read_count,fastq_ftp,submitted_ftp,sample_title&format=tsv&download=true&limit=0", 
              destfile = f2)


# the logic is to match required metadata files and only download
# the 16S rRNA files neccessary.  

metadata <- read_tsv(file = file.path(f1, "agp_metadata_20220304.txt"))
manifest <- read_table(file = f2)

subset_metadata <- metadata %>% select(sample_name, country_of_birth, 
                    age_years, description, host_body_site, host_taxid, 
                    diabetes, sex, ibd, ibs, ibd_diagnosis, ibd_diagnosis_refined) %>% 
    filter(host_body_site == "UBERON:feces") %>% 
    filter(description %in% c("American Gut Project Stool Sample", "American Gut Project Stool sample"), 
           host_taxid == "9606") %>% 
    select(-c(host_body_site, description, host_taxid))
joint_ids <- subset_metadata$sample_name
subset_manifest <- manifest %>% filter(library_strategy == "AMPLICON") %>% 
    filter(read_count >= 1000) %>% 
    filter(sample_title %in% subset_ids) %>% 
    select(sample_title, read_count, fastq_ftp, submitted_ftp) %>% 
    rename("sample_name" = "sample_title")

joint_metadata <- inner_join(subset_metadata, subset_manifest)

saveRDS(joint_metadata, file = "metadata/agp_joint_mtd.rds")



