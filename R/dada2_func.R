#' Filter samples based on metadata for gevers et al
#' @return file path for input data
filter_samples <- function(fastq_path, metadata_path) {
    files <- list.files(fastq_path)
    metadata <- read.table(file = metadata_path, sep = "\t", header = TRUE)
    sample_names <- map_chr(str_split(files, pattern = ".fastq.gz"), ~ {
        .x[1]
    })
    metadata <- metadata |>
        filter(biopsy_location == "Terminal ileum") |>
        select(sample_name, diagnosis) |>
        filter(diagnosis %in% c("no", "CD"))
    common <- intersect(sample_names, sample_metadata$sample_name)
    # get new path
    path_split <- strsplit(path, split = "/")[[1]]
    path_split <- path_split[1:(length(path_split) - 1)]
    path_split <- c(path_split, "input_data")
    input_path <- do.call(file.path, as.list(path_split))
    dir.create(input_path, showWarnings = FALSE, recursive = TRUE)
    # copy files over
    for (i in list.files(fastq_path, full.names = TRUE)) {
        names <- strsplit(strsplit(i, "//")[[1]][2], ".fastq.gz", fixed = TRUE)[[1]]
        if (names %in% common) {
            file.copy(i, input_path)
        }
    }
    return(input_file)
}


#' retrieve sequences and filter and trim
#' Assumes that this is a single-end files 
filter_and_trim <- function(reads, base_path=NULL, pat = NULL, control=NULL) {
    args <- list(
        truncLen=145,
        maxN=0,
        maxEE=2, 
        truncQ=2, 
        compress=TRUE,
        multithread=TRUE, 
        rm.phix = TRUE
    )
    if (!is.null(control)){
        if (!is(control, "list")){
            stop("Requires control arguments to be list")
        }
        ref_names <- c("truncLen", "maxN", "maxEE", "truncQ", 
                       "compress", "multithread", "rm.phix")
        if (length(setdiff(names(control), ref_names)) >= 1){
            warnings
        }
        i_names <- intersect(names(control), names(args))
        for (i in i_names){
            args[[i]] <- control[[i]]
        }
    }
    
    sample_names <- sapply(strsplit(basename(reads), pat), `[`, 1)
    filt_path <- file.path(base_path, "filtered", paste0(sample_names, "_filt.fastq.gz"))
    
    args_ftr <- c(args, list(fwd = reads, filt = filt_path))
    
    out <- do.call(dada2::filterAndTrim, args_ftr)
    return(filt_path)
}

#' learn_errors require output from filter_and_trim function
learn_errors <- function(filt_path) {
    print(filt_path)
    err <- dada2::learnErrors(filt_path, multithread = TRUE, 
                              nbases = 1e8, randomize = TRUE)
    return(err)
}

#' Run dada2 on learned errors and return outputs with chimeric sequences removed
run_dada2 <- function(filt_path, err) {
    derep <- derepFastq(filt_path)
    dada_samples <- dada2::dada(derep, 
                                err = err, 
                                multithread = TRUE)
    sample_name <- sapply(strsplit(basename(filt_path), "_filt.fastq.gz"), `[`, 1)
    print(sample_name)
    out <- list(dada_samples)
    names(out) <- sample_name
    return(out)
}

#' @param dada_res Dada res is a list of large size 
remove_chimera <- function(dada_res){
    # preprocess dada_res
    sample_names <- map_chr(dada_res, ~names(.x))
    dada_res <- map(dada_res, pluck,1)
    names(dada_res) <- sample_names
    seqtab <- dada2::makeSequenceTable(dada_res)
    seqtab_nochim <- dada2::removeBimeraDenovo(seqtab, method = "consensus", 
                                               multithread = TRUE, verbose = TRUE)
    return(seqtab_nochim)
}


#' Assign taxonomy
assign_taxonomy <- function(seqtab) {
    dir.create("databases", showWarnings = FALSE)
    if (!file.exists(here("databases", "silva_species_assignment_v138.1.fa.gz"))) {
        download.file("https://zenodo.org/record/4587955/files/silva_species_assignment_v138.1.fa.gz?download=1",
            destfile = "databases/silva_species_assignment_v138.1.fa.gz"
        )
    }
    if (!file.exists(here("databases","silva_nr99_v138.1_train_set.fa.gz"))) {
        download.file("https://zenodo.org/record/4587955/files/silva_nr99_v138.1_train_set.fa.gz?download=1",
            destfile = "databases/silva_nr99_v138.1_train_set.fa.gz"
        )
    }
    taxa <- assignTaxonomy(seqtab, 
                           here("databases", "silva_nr99_v138.1_train_set.fa.gz"), 
                           multithread = TRUE)
    taxa <- addSpecies(taxa, here("databases", "silva_species_assignment_v138.1.fa.gz"))
    output <- list(seqtab = seqtab, taxa = taxa)
    return(output)
}
