library(phyloseq)
library(tidyverse)

path_abun <- read_tsv(file = "data/pathabundances_3.tsv")
tax_abun <- read_tsv(file = "data/taxonomic_profiles_3.tsv")
metadata <- read_csv(file = "metadata/hmp2_metadata.csv")

metadata |> filter(data_type == "metagenomics") |>
   group_by(diagnosis, week_num) |> tally() |> arrange(-n) |> 
  filter(diagnosis == "nonIBD")
