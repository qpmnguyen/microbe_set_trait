# Functions commonly used  
library(mia)
library(taxizedb)
library(curatedMetagenomicData)
library(DBI)
library(ggtree)


hmp2 <- curatedMetagenomicData(pattern = "HMP_2012.relative_abundance", dryrun = FALSE)[[1]]
taxa <- rowData(hmp2)
rownames(taxa) <- paste0("tax", 1:nrow(hmp2))
tax_tree <- toTree(data = taxa)
rownames(hmp2) <- paste0("tax", 1:nrow(hmp2))
hmp2_new <- changeTree(x = hmp2, rowTree = tax_tree, rowNodeLab = taxa[["Species"]])



#' @title Function to take a string and initiate a query 
#' @param string The string of interest, should be of form k__Kingdom|p_Phylum 
#' @return An NCBI id as string 
query_db <- function(string){
  
  
}


#' @title Translating taxa table names to NCBI ids 
#' @description Access relevant databases to assign NCBI ids to taxa
#' @param seq (TreeSummarizedExperiment type). The sequence file to manipulate  
#' @param metaphlan (logical). To indicate whether the physeq is metaphlan results  
translate_ncbi <- function(seq, metaphlan = FALSE){
  if (class(tax_table != "TreeSummarizedExperiment")){
    stop("tax_table needs to be of TreeSummarizedExperiment class from package of the same name (or use mia)")
  }
  
  if (metaphlan == TRUE){
    pmap()
    
  }
}

