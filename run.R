#' Master script to run all individual pipelines
library(targets)
library(optparse)
library(tarchetypes)

option_list <- list(
    make_option(c("-c", "--ncores"), type="integer", default=5,
                help="Number of cores",
                metavar="NCORES"),
    make_option(c("-a", "--analysis"), type = "character",
                help = "What type of analysis to perform",
                metavar="ANALYSIS"),
    make_option(c("-r", "--remove"), type = "logical", default=TRUE,
                help = "Restart pipeline from scratch", metavar = "REMOVE"),
    make_option(c("-p", "--parallel"), type = "logical", default = TRUE, 
                help = "Run pipeline in parallel", metavar="PARALLEL")
)

opt <- parse_args(OptionParser(option_list=option_list))

name <- opt$analysis

match.arg(name, c("dada2_agp", "db_preprocess", "coverage", "pred"))

if (opt$remove == TRUE){
    tar_destroy()
}

Sys.setenv(TAR_PROJECT = name)

if (opt$parallel == TRUE){
    tar_make_future(workers = opt$ncores)
} else {
    tar_make()
}

