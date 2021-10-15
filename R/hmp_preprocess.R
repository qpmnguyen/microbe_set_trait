library(phyloseq)
library(tidyverse)

path_abun <- read_tsv(file = "data/pathabundances_3.tsv")
tax_abun <- read_tsv(file = "data/taxonomic_profiles_3.tsv")
metadata <- read_csv(file = "metadata/hmp2_metadata.csv")

metadata <- metadata |> filter(data_type == "metagenomics") |>
  filter(week_num == 0) |> select(`External ID`, diagnosis) |> 
  rename("ids" = `External ID`)

path_abun <- path_abun |> rename("path" = "Feature\\Sample") 

str_detect(string = path_abun$path[1])

           