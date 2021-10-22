# Functions commonly used  
library(mia)
library(taxizedb)
library(curatedMetagenomicData)
library(stringr)
library(purrr)
library(DBI)
library(ggtree)
library(data.table)


hmp2 <- curatedMetagenomicData(pattern = "HMP_2012.relative_abundance", dryrun = FALSE)[[1]]
taxa <- rowData(hmp2)
rownames(taxa) <- paste0("tax", 1:nrow(hmp2))
tax_tree <- toTree(data = taxa)
rownames(hmp2) <- paste0("tax", 1:nrow(hmp2))
hmp2_new <- changeTree(x = hmp2, rowTree = tax_tree, rowNodeLab = taxa[["Species"]])

test <- rownames(hmp2)[1:10]

#' @title Function to take a string and initiate a query 
#' @param string The string of interest, should be of form k__Kingdom|p_Phylum 
#' @return An NCBI id as string 
query_metaphlan <- function(str_vec){
  split_list <- str_split(string = str_vec, pattern = "\\|")
  query_list <- map(split_list, ~{
    sub <- gsub(.x, pattern = "[a-z]__", replacement = "")
    if (length(sub) != 8){
      query <- c(sub, rep(NA, 8 - length(sub)))
    }
    return(query)
  })
  rank_names <- c("Kingdom", "Phylum", "Class", "Order", 
                  "Family", "Genus", "Species", "Strain")
  res <- map(query_list, ~{
    q_list <- as.data.table(t(.x))
    names(q_list) <- rank_names
    q_df <- query_table[J(q_list), on = rank_names]
    q_df[,ncbiID,]
  })
  res <- flatten_chr(res)
  return(res)
}

#' @title Query standard NCBI database for names  
#' @param str_vec Vector of names 
#' @description Query only the last rank of the name 
query_standard <- function(str_vec, rank){
  
}


#' @title Translating taxa table names to NCBI ids 
#' @description Access relevant databases to assign NCBI ids to taxa
#' @param seq The sequence file to manipulate  
#' @param metaphlan (logical). To indicate whether the physeq is metaphlan results  
translate_ncbi <- function(seq, metaphlan = FALSE){
  if (class(seq) == "phyloseq"){
    
  } else if (class(seq) == "TreeSummarizedExperiment"){
    
  } else {
    stop("seq has to be phyloseq or TreeSummarizedExperiment")
  }
}

