library(phyloseq)
library(tidyverse)
library(data.table)
library(dtplyr)
library(taxizedb)
library(jsonlite)

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
  # replace sample names by removing _profile 
  old <- colnames(tax_abun)
  new <- c("tax", flatten_chr(strsplit(old[-1], split = "_profile")))
  setnames(tax_abun, old, new)
  # get only bacteria and the species level information 
  tax_abun <- tax_abun[str_detect(string = tax, pattern = "k__Bacteria") & str_detect(string = tax, pattern = "s__"), , ]
  # transpose
  tax_abun <- transpose(tax_abun, keep.names = "ids", make.names = "tax")
  # match to metadata
  metadata <- proc_metadata()
  tax_abun <- tax_abun[ids %in% metadata[,ids,],,]
  
  
  # reformat tax names
  taxtab <- matrix(0, nrow = ncol(tax_abun) - 1, ncol = 7)
  colnames(taxtab) <- c("Kingdom", "Phylum", "Class", "Order", 
                        "Family", "Genus", "Species")
  
  for (i in seq_len(ncol(tax_abun) - 1)){
    string_name <- colnames(tax_abun)[i+1]
    print(string_name)
    proc_names <- map_chr(flatten_chr(str_split(string_name, pattern = "\\|")), ~{
      gsub(x = .x, pattern = "[a-z]__", "")
    })
    taxtab[i,] <- proc_names
  }
  
}




