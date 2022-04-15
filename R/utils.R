# Functions commonly used
library(mia)
library(taxizedb)
library(curatedMetagenomicData)
library(stringr)
library(purrr)
library(DBI)
library(ggtree)
library(data.table)
library(phyloseq)
library(usethis)

#' @title Function to take a string and initiate a query
#' @param string The string of interest, should be of form k__Kingdom|p_Phylum
#' @return An NCBI id as string
query_metaphlan <- function(physeq) {
    if (class(physeq) == "phyloseq") {
        taxtab <- tax_table(physeq)
    } else {
        taxtab <- rowData(physeq)
    }

    query_list <- map(seq_len(nrow(taxtab)), function(x) {
        q_names <- as.vector(unname(as.matrix(taxtab[x, ])))
        if (length(q_names) < 8) {
            q_names <- c(q_names, rep(NA_character_, 8 - length(q_names)))
        }
        return(q_names)
    })

    query_table <- as.data.table(readRDS(file = file.path("metadata", "mpa_marker.rds")))
    # following is some old code w/r/t parsing the Metaphlan2 strings directly instead
    # split_list <- str_split(string = str_vec, pattern = "\\|")
    # query_list <- map(split_list, ~{
    #  sub <- gsub(.x, pattern = "[a-z]__", replacement = "")
    #  if (length(sub) != 8){
    #    query <- c(sub, rep(NA, 8 - length(sub)))
    #  }
    #  return(query)
    # })
    rank_names <- c(
        "Kingdom", "Phylum", "Class", "Order",
        "Family", "Genus", "Species", "Strain"
    )
    res <- map(query_list, ~ {
        q_list <- as.data.table(t(.x))
        names(q_list) <- rank_names
        q_list <- as.data.table(map(q_list, ~ {
            gsub(
                x = .x,
                pattern = " ",
                replacement = "_"
            )
        }))

        q_df <- query_table[J(q_list), on = rank_names]
        q_df[, ncbiID, ]
    })
    res <- flatten_chr(res)
    return(res)
}

#' @title Query standard NCBI database for names
#' @param str_vec Vector of names
#' @description Query only the last rank of the name
query_standard <- function(physeq, t_rank) {
    # first extract taxtab
    if (class(physeq) == "phyloseq") {
        taxtab <- tax_table(physeq)
        taxtab <- S4Vectors::as.data.frame(taxtab)
    } else {
        taxtab <- rowData(physeq)
    }
    # get the relevant rank names from ncbiids
    ranks <- c("superkingdom", "phylum", "class", "order", "family", "genus", "species")
    rank_seq <- seq_len(which(ranks == tolower(t_rank)))
    ranks <- ranks[rank_seq]

    # map over all the rows of taxtab
    ncbi_ids <- map_chr(seq_len(nrow(taxtab)), function(x) {
        ids <- taxizedb::name2taxid(taxtab[x, t_rank], out_type = "summary") |> dplyr::pull(id)
        # if ambiguous compare the entire name series
        if (length(ids) >= 2) {
            full_rank <- taxizedb::classification(ids)
            # compare the name series
            # so far no way to solve for heterotypic synonyms except to rely on the initial NCBI query
            rank_test <- map_lgl(full_rank, function(y) {
                ref_name <- as.vector(as.matrix(unname(taxtab[x, ])))
                new_name <- y |>
                    dplyr::filter(rank %in% ranks) |>
                    dplyr::pull(name)
                all(new_name == ref_name)
            })
            # if none of the entire series match, then return everything
            # assuming ncbi calls have solved the naming issues
            if (sum(rank_test) == 0) {
                return(ids)
            } else {
                return(ids[rank_test])
            }
            # return NAs if nothing returns
        } else if (length(ids) == 0) {
            return(NA_character_)
        } else {
            return(ids)
        }
    })
    return(ncbi_ids)
}


#' @title Translating taxa table names to NCBI ids
#' @description Access relevant databases to assign NCBI ids to taxa
#' @param seq The sequence file to manipulate
#' @param metaphlan (logical). To indicate whether the physeq is metaphlan results
translate_ncbi <- function(seq, t_rank = "genus", metaphlan = FALSE) {
    if (metaphlan == TRUE) {
        # the marker database only has absolute information for species
        if (t_rank != "Species") {
            stop("Rank has to be species to query to metaphlan")
        }
        queried_names <- query_metaphlan(seq)
    } else {
        queried_names <- query_standard(seq, t_rank = t_rank)
    }
    print(queried_names)
    if (class(seq) == "phyloseq") {
        taxa_names(seq) <- queried_names
    } else if (class(seq) == "TreeSummarizedExperiment") {
        rownames(seq) <- queried_names
    } else {
        stop("seq has to be phyloseq or TreeSummarizedExperiment")
    }
    return(seq)
}

#' @param name Name of the package 
#' @param version (logical). Indicate whether versions are kept together 
record_package <- function(name, version=TRUE){
    if (version){
        version <- packageVersion(name)
    } else {
        version <- NULL
    }
    usethis::use_package(name, type = "Imports", min_version=version)
}


