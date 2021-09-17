library(tidyverse)

weissman <- read_csv(file = "data/weissman.csv")
madin_traits <- read_csv(file = "data/madin_condensed_traits_NCBI.csv")
madin_species <- read_csv(file = "data/madin_condensed_species_NCBI.csv")

metadata <- colnames(weissman)[1:19]

w_db <- weissman |> 
  select(c(all_of(metadata),starts_with(c("Enzyme", "Volatile", "Substrate")))) |>
  pivot_longer(starts_with(c("Enzyme", "Volatile", "Substrate")), 
               names_to = "names", values_to = "values") |> 
  separate(names, c("class", "trait"), sep = "\\.\\.")


metadata <- colnames(madin_traits)[1:11]
mt_db <- madin_traits |> 
  select(c(all_of(metadata), c("metabolism", "pathways", "carbon_substrates")))

metadata <- colnames(madin_species)[1:9]
ms_db <- madin_species |> 
  select(c(all_of(metadata), c("metabolism", "pathways", "carbon_substrates")))

ms_db |> pull(carbon_substrates) |> unique()

substrates <- ms_db |> pull(carbon_substrates) |> unique()

substrates <- sapply(substrates, function(x) str_split(x, pattern = ","))

substrates <- flatten_chr(substrates) |> unique()


pathways <- ms_db |> pull(pathways) |> unique()

pathways <- sapply(pathways, function(x) str_split(x, pattern = ","))

pathways <- flatten_chr(pathways) |> unique()
