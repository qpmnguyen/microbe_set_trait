library(phyloseq)
library(tidyverse)
library(data.table)
library(dtplyr)

proc_metadata <- function(){
  metadata <- fread("metadata/hmp2_metadata.csv")
  metadata <- metadata[data_type == "metagenomics" & week_num == 0,.(`External ID`, diagnosis)]
  setnames(metadata, "External ID", "ids")
  return(metadata)
}

proc_pathabun <- function(){
  path_abun <- fread("data/pathabundances_3.tsv.gz")
  old <- colnames(path_abun)
  new <- c("path", flatten_chr(strsplit(old[-1], split = "_pathabundance_cpm")))
  setnames(path_abun, old, new)
  path_abun <- path_abun[!str_detect(string = path, pattern = "\\|"),,]
  path_abun <- transpose(path_abun, keep.names = "ids", make.names = "path")
  metadata <- proc_metadata()
  path_abun <- path_abun[ids %in% metadata[,ids,], , ]
  return(path_abun)
}



proc_taxabun <- function(){
  tax_abun <- fread("data/taxonomic_profiles_3.tsv.gz")
  old <- colnames(path_abun)
  new <- c("path", flatten_chr(strsplit(old[-1], split = "_pathabundance_cpm")))
  colnames(tax_abun)[1:10]
}