library(tidyverse)
library(RSQLite)
library(DBI)
library(jsonlite)

filepath <- file.path("large_files", "mpa_v30_CHOCOPhlAn_201901_marker_info.txt")

#' @title Processing metaphlan3 marker information using json
#' @param string The string to the processed as a json file
#' @return A data frame containing the individual taxonomic levels as specified
process_json <- function(string) {
    string <- gsub(string, pattern = "'", replacement = "\"")
    string_l <- jsonlite::fromJSON(txt = string)
    taxon <- string_l$taxon
    split <- str_split(taxon, pattern = "\\|")[[1]]
    df_taxon <- data.frame(t(rep(NA, 8)))
    colnames(df_taxon) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain")

    p_split <- gsub(split, pattern = "[a-z]__", replacement = "")
    if (length(p_split) != ncol(df_taxon)) {
        p_split <- c(p_split, rep(NA, ncol(df_taxon) - length(p_split)))
    }
    df_taxon[1, ] <- p_split
    return(df_taxon)
}


# code for function to read in a file streaming in R
# https://stackoverflow.com/questions/12626637/read-a-text-file-in-r-line-by-line
#' @title Create metaphlan marker db
#' @description This function takes in the mpa marker information file from metaphlan3
#'  and create an equivalent sql data base to search through.
#' @param filepath Path to the marker information txt file
#' @param outpath Path to the resulting sqlite database.
create_db <- function(filepath, outpath = "metadata/mpa_marker.db") {
    dbcon <- dbConnect(RSQLite::SQLite(), outpath)
    # if there is no table, create table
    if (length(dbListTables(dbcon)) == 0) {
        fields <- c("INTEGER", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT", "TEXT")
        names(fields) <- c("ncbiID", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain")
        dbCreateTable(dbcon, name = "taxonomy", fields = fields)
    }

    # stream through file
    filecon <- file(filepath, "r")
    count <- 1
    while (TRUE) {
        line <- readLines(filecon, n = 1)
        if (length(line) == 0) {
            break
        }
        proc_line <- str_split(line, pattern = "\t")[[1]]
        ncbiid <- as.integer(str_split(proc_line[1], "__")[[1]][1])
        # search to see if the queried taxa is already in the database
        search_id <- dbGetQuery(dbcon, "SELECT * FROM taxonomy where ncbiid = ?")
        dbBind(search_id, list(ncbiid))
        queried_id <- dbFetch(search_id)
        if (nrow(queried_id) == 0) {
            dbClearResult(search_id)
            insert_line <- process_json(proc_line[2])
            insert_line <- cbind(ncbiid, insert_line)
            dbAppendTable(conn = dbcon, name = "taxonomy", value = insert_line)
            count <- count + 1
            if (count %% 500 == 0) {
                print(count)
            }
        } else {
            next
        }
    }
    close(filecon)
    dbDisconnect(dbcon)
}
