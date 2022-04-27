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
    c_ranks <- c(
        "superkingdom", "phylum", "class",
        "order", "family", "genus"
    )
    query <- classification(spec_id)[[1]]
    if (!is.data.frame(query)){
        output <- list(
            id = NA_character_,
            full_name = NA_character_
        )
    } else {
        output <- list(
            id = query %>% filter(rank == "genus") %>% pull(id),
            full_name = query %>% filter(rank %in% c_ranks) %>% 
                pull(name) %>% str_c(collapse = "_")
        )
        query <- query %>% filter(rank == "genus")
    }
    return(output)
}