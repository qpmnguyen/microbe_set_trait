library(tidyverse)

weissman <- read_csv(file = "data/weissman.csv")

head(weissman)

metadata <- colnames(weissman)[1:19]

w_db <- weissman |> 
  select(c(all_of(metadata),starts_with(c("Enzyme", "Volatile", "Substrate")))) |>
  pivot_longer(starts_with(c("Enzyme", "Volatile", "Substrate")), 
               names_to = "names", values_to = "values") |> 
  separate(names, c("class", "trait"), sep = "\\.\\.")

