library(tidyverse)
library(phyloseq)
library(glue)
library(ape)
library(Biostrings)

#' Load the gevers et al data set and reformat 
load_and_format <- function(file_path = "data/gevers_dada2.rds"){
  data <- readRDS(file = file_path)
  otu_tab <- data$seqtab_nochim
  tax_tab <- data$taxa
  
  # getting sequences 
  seq <- colnames(otu_tab)
  names(seq) <- paste0("ASV", 1:length(seq))
  seq <- Biostrings::DNAStringSet(seq)
  colnames(otu_tab) <- names(seq)
  rownames(tax_tab) <- names(seq)
  
  # getting metadata 
  metadata <- read_tsv(file = "metadata/gevers_metadata.txt")
  samp_names_fmt <- map_chr(str_split(rownames(otu_tab), 
                                      pattern = ".fastq.gz"), ~{ .x[1]})
  rownames(otu_tab) <- samp_names_fmt
  metadata <- metadata |> filter(sample_name %in% samp_names_fmt) |> 
    select(sample_name, age, diagnosis) |> column_to_rownames(var = "sample_name")
  
  # stitch together using phyloseq
  physeq <- phyloseq(
    otu_table(otu_tab, taxa_are_rows = FALSE),
    sample_data(metadata),
    tax_table(tax_tab),
    seq
  )
  
  # let's do some filtering!  
  #physeq <- physeq |> 
  #  filter_taxa(function(x) sum(x > 0) > (0.1 * length(x)), TRUE)
  
  return(physeq)
}

export_picrust <- function(physeq){
  Biostrings::writeXStringSet(refseq(physeq), filepath = "data/gevers_seq.fa")
  otu <- t(as(otu_table(physeq), "matrix"))
  otu_biom <- biomformat::make_biom(data = otu)
  biomformat::write_biom(otu_biom, "data/gevers.biom")
}


perform_picrust <- function(nthreads = 2){
  if (!file.exists("picrust2.nf")){
    stop("Require nextflow command")
  } 
  if (!"picrust_2.4.1" %in% reticulate::conda_list()$name){
    stop("Please create environment named picrust_2.4.1 from picrust_env.yml file")
  }
  cmd <- glue::glue("source ~/.bashrc && conda activate picrust_2.4.1 && nextflow run picrust2.nf --threads {nthreads}", nthreads = nthreads)
  print(cmd)
  system(cmd)
}