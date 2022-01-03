# scripts to pre-process and match the metacyc and
# taxonomic databases
library(data.table)
library(magrittr)
library(phyloseq)
library(BiocSet)
library(taxizedb)
library(glue)
library(taxadb)
library(purrr)
library(DBI)
library(stringr)
library(dplyr)
library(dbplyr)
library(tidyr)

if (length(tdb_cache$list()) == 0) {
    taxizedb::db_download_ncbi(verbose = FALSE, overwrite = FALSE)
}

# PROCESS THE TRAIT DATABASE ####   

#' Read in trait database as a data.table  
read_trait <- function(){
    trait_db <- fread(file = "data/condensed_species_NCBI.txt")
    return(trait_db)
}

#' subset database by trait
#' @param traitdb The trait database as a data.table 
#' @param trait A string indicating the desired trait to be extracted 
get_subset <- function(trait_db, trait) {
    trait <- match.arg(trait, c("pathways", "carbon_substrates", 
                                "sporulation", 
                                "gram_stain", "cell_shape", "range_tmp", 
                                "range_salinity", 
                                "motility", "metabolism"))
    metadata <- c(colnames(trait_db)[1:8])
    columns <- c(metadata, trait)
    subset_db <- trait_db[,..columns] %>% .[superkingdom == "Bacteria",,] %>% 
        .[!is.na(i), ,env = list(i = trait)]
    return(subset_db)
}

#' @title Retrieve ncbiids based on unique traits 
#' @param subset_db The trait database ideally restricted to 
#'     the trait of interest 
#' @param trait String indicating the trait of interest 
get_traits <- function(subset_db, trait){
    trait <- match.arg(trait, c("pathways", "carbon_substrates", 
                                "sporulation", 
                                "gram_stain", "cell_shape", "range_tmp", 
                                "range_salinity", 
                                "motility", "metabolism"))
    idx <- ncol(subset_db)
    if (trait %in% c("pathways", "carbon_substrates")){
        subset_db <- subset_db[,j := map(str_split(j, ","), str_trim), 
                               env = list(j = trait)]
        trait_list <- subset_db[,j, env = list(j = trait)] %>% 
            flatten_chr() %>% unique()
    } else {
        trait_list <- subset_db[,j, env = list(j = trait)] %>% unique()
    }

    
    # go through each trait in the trait list, retrieve a list of NCBI ids  
    # and then return everything  
    ncbiid_list <- map(trait_list, function(x){
        subset_db[map_lgl(i, ~{x %in% .x}),,env = list(i = trait)] %>% 
            pull(species_tax_id)
    })
    names(ncbiid_list) <- trait_list
    return(ncbiid_list)
}


#' @title Getting the biocSets based on whether aggregation should be performed 
#' @param ncbiid_list The list of ncbiids 
#' @param g_agg Logical indicating whether aggregating to genus level 
get_sets <- function(ncbiid_list, trait_db, trait, g_agg=TRUE) {
    if (g_agg == FALSE){
        ncbiid_list <- map(ncbiid_list, as.character)
        sets <- BiocSet::BiocSet(ncbiid_list)
        sizes <- es_elementset(sets) %>% group_by(set) %>% count()
        sets <- sets %>% mutate_set(type = trait, size = sizes$n)
    } else {
        

        
    }
    
    
    return(0)
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
        # hypergeometric distribution test for underrepresentation 
        # P(X <= x) 
        cum_prob <- phyper(
            q = db_ngenus,
            m = ncbi_ngenus,
            n = ncbi_nspec - ncbi_ngenus,
            k = db_nspec,
            lower.tail = TRUE
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



#' @title Function to perform testing for each set of traits within a trait
#'     type 
#' @param id_vec The identifier ncbiids in vector format (NOT A LIST)
#' @param trait_db The trait database as data.table format
get_filtered_genus <- function(id_vec, trait_db){
    c_ranks <- c(
        "superkingdom", "phylum", "class",
        "order", "family", "genus"
    )
    # first let's get all the counts for the species within a genera 
    # with the trait  
    tr_counts <- trait_db[species_tax_id %in% id_vec, .(.N), by = c_ranks]
    
    # second let's get all the counts for the genera available in the database 
    tot_counts <- trait_db[,.(.N), by = c_ranks]
    
    # merge them together and calculate proportions  
    merged_counts <- tot_counts[tr_counts, by = c_ranks][,prop := i.N/N]
    
    
    
    
}


# PROCESS NAMES INTO NCBIIDS ####

# call ncbi ids using taxizedb function 
# Use the clean_names function from taxadb to clean names 
#' @param vec A vector of names represented 
#'     by ranks from superkingdom to genus
call_genus_id <- function(vec) {
    c_ranks <- c(
        "superkingdom", "phylum", "class",
        "order", "family", "genus"
    )
    # clean names and extract genus
    vec <- as.matrix(vec)
    vec <- map_chr(vec, ~taxadb::clean_names(.x))
    genus <- vec[length(vec)]
    
    # retrieve ncbiids
    id_names <- name2taxid(genus, out_type = "summary")
    
    # if there are ambiguous idenfiers
    if (nrow(id_names) >= 2){
        message("Ambiguous identifiers")
        # compare reference with new names and return a logical vector
        # equal to the length of each match  
        ref_name <- vec
        check_names <- map_chr(id_names %>% pull(id), ~{
            query <- classification(.x, db = "ncbi")[[1]]
            if (is.data.frame(query) == FALSE){
                return(NA_character_)
            } else {
                new_name <- query %>% 
                    dplyr::filter(rank %in% c_ranks) %>% 
                    dplyr::pull(name)
                matching <- ref_name == new_name 
                diff <- length(matching) - length(which(matching))
                # handling cases where there are mismatches
                if (diff <= 1) {
                    return(TRUE)
                } else {
                    return(FALSE)
                }
            }
        })
        # filter out all false 
        id_names <- id_names[which(check_names == TRUE & !is.na(check_names))]
        # if there is still length larger than 1, return the top option 
        if (nrow(id_names) >= 2){
            id_names <- id_names$id[1]
        } else {
            id_names <- id_names$id
        }
    } else {
        id_names <- id_names$id
    }
    return(id_names)
}


#' @title Function to use the species name to back-track genera names  
#' @param spec_id Species NCBI identifier  
genus_from_species <- function(spec_id){
    query <- classification(spec_id)[[1]]
    if (!is.data.frame(query)){
        output <- list(
            id = NA_character_,
            full_name = NA_character_
        )
    } else {
        output <- list(
            id = query %>% filter(rank == "genus") %>% pull(id)
            full_name = query %>% filter(rank == )
        )
        query <- query %>% filter(rank == "genus")
    }
    return(output)
}





# PROCESS NCBIID RESULTS ####  
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