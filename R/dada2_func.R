library(dada2)
library(ShortRead)
library(tidyverse)
library(rlist)

#' retrieve sequences and filter and trim
filter_and_trim <- function(path){
  fnFs <- sort(list.files(path, pattern="fastq.gz", full.names = TRUE))
  sample_names <- sapply(strsplit(basename(fnFs), ".fastq.gz"), `[`, 1)
  filt_path <- file.path(path, "filtered")

  out <- dada2::filterAndTrim(fnFs, filt_path, 
                       truncLen = 165, maxN = 0, maxEE = 2, truncQ = 2, rm.phix = TRUE, 
                       compress = TRUE, multithread = TRUE)
  output <- list(
    path = path,
    sample_names = sample_names,
    out = out, 
    filt_path = filt_path
  )
  return(output)
}

#' learn_errors require output from filter_and_trim function 
learn_errors <- function(output){
  if (!setequal(names(output), c("path", "sample_names", "out", "filt_path"))){
    stop("Output does not have all the required elements!")
  }
  if (length(list.files(output$filt_path, pattern = "fastq.gz")) != length(list.files(output$path, pattern = "fastq.gz"))){
    stop("No files exist!")
  }
  err <- dada2::learnErrors(output$filt_path, multithread = TRUE)
  output <- rlist::list.append(output, err = err)
  return(output)
}

#' Run dada2 on learned errors and return outputs with chimeric sequences removed 
run_dada2 <- function(output){
  if (!setequal(names(output), c("err", "path", "sample_names", "out", "filt_path"))){
    stop("Output does not have all the required elements!")
  }
  dada_samples <- dada2::dada(output$filt_path, err = output$err, multithread = TRUE)
  seqtab <- dada2::makeSequenceTable(dada_samples)
  seqtab_nochim <- dada2::removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
  
  # track reads through the pipeline 
  getN <- function(x) sum(getUniques(x))
  track <- cbind(output$out, sapply(dada_samples, getN), rowSums(seqtab), rowSums(seqtab_nochim))
  # If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
  colnames(track) <- c("input", "filtered", "denoised", "w_chimeria", "nonchim")
  rownames(track) <- output$sample_names
  output <- rlist::list.append(output, 
                               track = track,
                               seqtab_nochim = seqtab_nochim)
  output <- rlist::list.remove(output, "out")
}

#' Assign taxonomy  
assign_taxonomy <- function(output){
  dir.create("databases", showWarnings = FALSE)
  if (!setequal(names(output), c("err", "path", "sample_names", "filt_path", "track", "seqtab_nochim"))){
    stop("Output does not have all the required elements!")
  }
  if (!file.exists("databases/silva_species_assignment_v138.1.fa.gz")){
    download.file("https://zenodo.org/record/4587955/files/silva_species_assignment_v138.1.fa.gz?download=1", 
                  destfile = "databases/silva_species_assignment_v138.1.fa.gz")
  }
  if (!file.exists("databases/silva_nr99_v138.1_train_set.fa.gz")){
    download.file("https://zenodo.org/record/4587955/files/silva_nr99_v138.1_train_set.fa.gz?download=1", 
                  destfile = "databases/silva_nr99_v138.1_train_set.fa.gz")
  }
  taxa <- assignTaxonomy(output$seqtab_nochim, "databases/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)
  taxa <- addSpecies(taxa, "databases/silva_species_assignment_v138.1.fa.gz")
  output <- rlist::list.append(output, taxa = taxa)
  return(output)
}





