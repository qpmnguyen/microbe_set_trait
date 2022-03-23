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
filter_and_trim <- function(path, pat="_001.fastq", control=NULL, limit=FALSE) {
    args <- list(
        truncLen=145,
        maxN=0,
        maxEE=2, 
        truncQ=11, 
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
    # first, extract reads
    reads <- sort(list.files(file.path(path, "raw"), pattern = pat, full.names = TRUE))
    
    if (limit) {
        reads <- reads[1:10]
    }
    
    sample_names <- sapply(strsplit(basename(reads), pat), `[`, 1)
    filt_path <- file.path(path, "filtered", paste0(sample_names, "_filt.fastq.gz"))
    
    args_ftr <- c(args, list(fwd = reads, filt = filt_path))
    
    out <- do.call(dada2::filterAndTrim, args_ftr)
    output <- list(
        path = path,
        sample_names = sample_names,
        out = out,
        filt_path = filt_path
    )
    return(output)
}

#' learn_errors require output from filter_and_trim function
learn_errors <- function(output) {
    if (!setequal(names(output), c("path", "sample_names", "out", "filt_path"))) {
        stop("Output does not have all the required elements!")
    }
    err <- dada2::learnErrors(output$filt_path, multithread = TRUE)
    output <- c(output, list(err = err))
    return(output)
}

#' Run dada2 on learned errors and return outputs with chimeric sequences removed
run_dada2 <- function(output) {
    if (!setequal(names(output), c("err", "path", "sample_names", "out", "filt_path"))) {
        stop("Output does not have all the required elements!")
    }
    dada_samples <- dada2::dada(output$filt_path, 
                                err = output$err, 
                                multithread = TRUE)
    seqtab <- dada2::makeSequenceTable(dada_samples)
    seqtab_nochim <- dada2::removeBimeraDenovo(seqtab, 
                                               method = "consensus", 
                                               multithread = TRUE, verbose = TRUE)

    # track reads through the pipeline
    getN <- function(x) sum(getUniques(x))
    track <- cbind(output$out, sapply(dada_samples, getN), rowSums(seqtab), rowSums(seqtab_nochim))
    # If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
    colnames(track) <- c("input", "filtered", "denoised", "w_chimeria", "nonchim")
    rownames(track) <- output$sample_names
    output <- rlist::list.append(output,
        track = track,
        seqtab_nochim = seqtab_nochim
    )
    output <- rlist::list.remove(output, "out")
}

#' Assign taxonomy
assign_taxonomy <- function(output) {
    dir.create("databases", showWarnings = FALSE)
    if (!setequal(names(output), c("err", "path", "sample_names", "filt_path", "track", "seqtab_nochim"))) {
        stop("Output does not have all the required elements!")
    }
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
    taxa <- assignTaxonomy(output$seqtab_nochim, 
                           here("databases", "silva_nr99_v138.1_train_set.fa.gz"), multithread = TRUE)
    taxa <- addSpecies(taxa, here("databases", "silva_species_assignment_v138.1.fa.gz"))
    output <- rlist::list.append(output, taxa = taxa)
    return(output)
}
