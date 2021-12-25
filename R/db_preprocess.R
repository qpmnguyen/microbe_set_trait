# scripts to pre-process and match the metacyc and
# taxonomic databases
library(tidyverse)
library(phyloseq)
library(BiocSet)
library(taxizedb)
library(reticulate)
library(glue)
library(dbplyr)

if (length(tdb_cache$list()) == 0) {
    taxizedb::db_download_ncbi(verbose = FALSE, overwrite = FALSE)
}

# call ncbi ids
#' @param str_ids A list of vectors of names from superkingdom to genus level
#' @param t_rank Unused argument for future control of which ranks to call on the ids
call_id <- function(str_ids, t_rank = "genus") {
    # get genus identifiers
    ids <- map_chr(str_ids, ~ {
        .x[length(.x)]
    })
    # map over genus names
    ncbi <- map(ids, ~ {
        if (is.na(.x)) {
            return(NA_character_)
        } else {
            res <- name2taxid(.x, out_type = "summary")
            # remove anything that's not bacteria
            ids <- map_chr(res$id, ~ {
                classif <- classification(.x)
                superkingdom <- classif[[1]] |>
                    filter(rank == "superkingdom") |>
                    pull(name)
                if (superkingdom == "Bacteria") {
                    return(.x)
                } else {
                    return(NA_character_)
                }
            })
            ids <- as.vector(na.omit(ids))
            # if after removing bacteria, there are still ids remain
            if (length(ids) >= 2) {
                # get the relevant ranks for the classification
                ranks <- c(
                    "superkingdom", "phylum", "class",
                    "order", "family", "genus", "species"
                )
                rank_seq <- seq_len(which(ranks == tolower(t_rank)))
                ranks <- ranks[rank_seq]
                full_rank <- classification(ids)
                # test for equivalence for each full rank name
                rank_test <- map_lgl(full_rank, function(y) {
                    ref_name <- str_ids
                    new_name <- y |>
                        dplyr::filter(rank %in% ranks) |>
                        dplyr::pull(name)
                    all(ref_name == new_name)
                })

                # if all rank_test is FALSE, return everything
                if (sum(rank_test) == 0) {
                    message("Ambiguous identifiers")
                    return(ids)
                } else {
                    return(ids[rank_test])
                }
            } else {
                return(ids)
            }
        }
    })
    return(ncbi)
}


#' @param trait_vec Vector of unique traits
#' @param db Database where species/taxa as rows and traits as columns.
#' @param filt_term Term to filter database on
#' @param id_col Name of id column indicating NCBI identifiers
#' @param genus_agg What are the direct levels
#' @param threshold Threshold for when a trait is considered
retr_sets <- function(db, trait_vec, filt_term,
                      id_col, genus_agg = FALSE, threshold = 0.95) {
    # special processing if filtered term is carbon_substrates, pathways
    db_nspec <- nrow(db)
    set_list <- map(trait_vec, ~ {
        if (filt_term %in% c("carbon_substrates", "pathways")) {
            db_inter <- db |>
                rowwise() |>
                mutate(t_bool = str_detect(
                    string = !!sym(filt_term),
                    pattern = .x
                )) |>
                ungroup()
        } else {
            db_inter <- db |>
                mutate(t_bool = ifelse(!!sym(filt_term) == .x, TRUE, FALSE))
        }
        if (genus_agg == TRUE) {
            c_ranks <- c(
                "superkingdom", "phylum", "class",
                "order", "family", "genus"
            )

            # first, filter out genera without trait at at least 0.95 of species
            print(.x)
            db_agg <- db_inter |>
                group_by(across(all_of(c_ranks))) |>
                summarise(t_count = sum(t_bool, na.rm = TRUE), n = n()) |>
                filter(t_count / n >= threshold)

            full_name <- pmap(db_agg |> select(-c(t_count, n)), ~ {
                unname(c(...))
            })
            # retreive ncbiids
            ncbi_ids <- call_id(full_name)

            # get p_values for filter
            p_values <- imap(ncbi_ids, ~ {
                test_genus(.x, db_nspec = db_nspec, db_ngenus = db_agg$n[.y])
            })
            sig_ids <- map(p_values, proc_gtest)
            flatten_chr(sig_ids)
        } else {
            db_inter |>
                filter(t_bool == TRUE) |>
                pull(!!sym(id_col)) |>
                as.character()
        }
    })
    names(set_list) <- trait_vec
    sets <- BiocSet(set_list)
    # remove NAs if there are any
    sets <- sets |> filter_element(!is.na(element))
    sets <- sets |> mutate_set(type = filt_term)
    return(sets)
}

#' This function returns a processed trait_set
process_db <- function() {
    trait_db <- read.table(
        file = "data/condensed_species_NCBI.txt",
        sep = ",", header = TRUE
    )
    metadata <- colnames(trait_db)[1:8]
    traits <- c("metabolism", "pathways", "carbon_substrates", "sporulation")
    # restrict databases to certain types of traits
    trait_db <- trait_db |>
        as_tibble() |>
        select(all_of(metadata), all_of(traits)) |>
        dplyr::filter(superkingdom == "Bacteria")

    return(trait_db)
}

#' This function takes the trait database and construct
get_trait_list <- function(trait_db) {
    traits <- c("metabolism", "pathways", "carbon_substrates", "sporulation")
    trait_list <- map(traits, ~ {
        t_vec <- unique(trait_db |> pull(!!sym(.x)))
        if (.x %in% c("pathways", "carbon_substrates")) {
            t_list <- map(t_vec, ~ str_split(.x, ",")) |>
                flatten() |>
                flatten_chr()
            t_vec <- map_chr(t_list, ~ str_trim(.x)) |>
                unique() |>
                na.omit() |>
                as.vector()
        } else {
            t_vec <- na.omit(t_vec)
        }
        return(t_vec)
    })
    return(trait_list)
}

#' @title Generate sets
#' @param trait_list List of traits prepared
#' @param trait_db The database prepared
#' @param agg Whether we're aggregating to Genus level
#' @param threshold Step 2 of aggregation threshold for a trait to be assigned
create_sets <- function(trait_list, trait_db, agg = FALSE, threshold = 0.95) {
    traits <- c("metabolism", "pathways", "carbon_substrates", "sporulation")
    trait_sets <- imap(trait_list, ~ {
        retr_sets(
            db = trait_db, trait_vec = .x, filt_term = traits[.y],
            id_col = "species_tax_id", genus_agg = agg, threshold = threshold
        )
    })
}

# do not run if trying to perform this reproducibly
# this is due to the fact that it requires an external software
# Please visit pathwaytools page for more information.
process_metacyc <- function() {
    # please make sure pathway tools is running in the background
    # The following commands did not work because it forces
    # the current R console to be in the pathway-tools program
    # with no easy way to exit
    # pthway_tools_exec <- "~/pathway-tools"
    # system(glue("{dir}/pathway-tools -lisp -python-local-only",
    #            dir = pthway_tools_exec))
    Sys.setenv(RETICULATE_PYTHON = "~/miniconda3/envs/microbe_set_trait/bin/python3")
    reticulate::source_python("python/metacyc.py")
    # get_hierarchy and get_instances are python functions
    subclass <- get_hierarchy()
    inst <- map(subclass, ~ {
        get_instances(.x)
    })
    names(inst) <- subclass
    return(inst)
}

#' n_spec is not part of the function since it is an expensive operation that
#' should be done only once
#' @param genus_id NCBI identifier for genus as string or numeric values. If
#'     string, convert to numeric
#' @param db_nspec Total number of species in the trait database, considered to be
#'     equivalent to the paramter k in a hypergeometric distribution
#' @param db_ngenus Total number of species assigned to the genus, considered to be
#'     equivalent to the parameter q in a hypergeometric distribution
#' @param ncbi_nspec Total number of species in NCBI. Considered to be m + n parameter (or nn)
#'     in a hypergeometric distribution
test_genus <- function(genus_id, db_nspec, db_ngenus) {
    if (!is.vector(genus_id)) {
        genus_id <- as.vector(genus_id)
    }
    genus_id <- map_dbl(genus_id, as.numeric)
    con <- DBI::dbConnect(RSQLite::SQLite(), taxizedb::tdb_cache$list()[1])
    ncbi_nspec <- tbl(con, sql("SELECT COUNT(*) FROM nodes WHERE rank = 'species'")) |>
        dplyr::pull()
    p_vals <- map_dbl(genus_id, function(x) {
        # get the table where the rank is species and the parent tax id is the genus
        # id that we've been talking about
        tbl <- tbl(con, sql("SELECT tax_id, parent_tax_id, rank FROM nodes"))
        # the total number of species assigned to that genus in ncbi database
        ncbi_ngenus <- tbl |>
            filter(rank == "species", parent_tax_id == x) |>
            tally() |>
            pull(n)
        # hypergeometric distribution
        cum_prob <- phyper(
            q = db_ngenus - 1,
            m = ncbi_ngenus,
            n = ncbi_nspec - ncbi_ngenus,
            k = db_nspec,
            lower.tail = FALSE
        )
        cum_prob
    })
    DBI::dbDisconnect(con)
    names(p_vals) <- genus_id
    return(p_vals)
}

# this function process a list of pvalues and returns NAs if
# things go wrong
proc_gtest <- function(pvals) {
    if (length(pvals) == 0) {
        return(NA_character_)
    }
    sig <- which(pvals < 0.05)
    if (length(sig) == 0) {
        return(NA_character_)
    } else {
        return(names(pvals)[sig])
    }
}
