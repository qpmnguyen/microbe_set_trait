library(tidyverse)
library(RSQLite)
library(DBI)
filepath <- "large_files/mpa_v30_CHOCOPhlAn_201901_marker_info.txt"
processFile = function(filepath) {
  count <- 0
  con = file(filepath, "r")
  while ( TRUE ) {
    line = readLines(con, n = 1)
    count <- count + 1
    if ( count == 10 ) {
      break
    }
    print(line)
  }
  close(con)
}

# code for function to read in a file streaming in R 
# https://stackoverflow.com/questions/12626637/read-a-text-file-in-r-line-by-line
create_db <- function(){
  dbcon <- dbConnect(RSQLite::SQLite(), "metadata/mpa_marker.db")
}