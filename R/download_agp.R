library(here)
library(tidyverse)

here::i_am("R/download_agp.R")
dir.create(here("raw","agp"))

manifest <- readRDS(file = here("metadata", "agp_joint_mtd.rds") 

purrr::map(seq_along(nrow(manifest)), ~{
    filename <- paste0(manifest$sample_name[.x], "_001.fastq")
    path <- here("raw", "agp", filename)
    download.file(url = manifest$fastq_ftp[.x], destfile = path)
    return(0)
i})
