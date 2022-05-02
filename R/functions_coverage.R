library(tidyverse)
library(phyloseq)
library(vegan)

#' This function counts up for each sample the number of sets having 
#' at least 2 present taxa represented  
count_traits <- function(physeq, q_set){
    eset <- q_set %>% es_elementset() %>% as_tibble() 
    tab <- otu_table(physeq) %>% as.matrix()
    tab[tab > 0] <- 1
    # obtain list of present otus 
    present <- apply(tab, 2, function(x){
        idx <- which(x >= 1)
    })
    present <- map(present, ~rownames(tab)[.x])
    present <- map_dbl(present, ~{
        eset %>% filter(element %in% .x) %>% 
            group_by(set) %>% dplyr::tally() %>% filter(n >= 2) %>% nrow()
    })
    return(present)
}

#' Count the diversity of taxa per sample per set 
count_tax <- function(physeq, q_set, mode){
    mode <- match.arg(mode, c("evenness", "richness"))
    # get all taxa members  
    p_tax <- q_set@element %>% pull(element)
    # filter out the physeq for only taxa that are part of the sets 
    red_physeq <- prune_taxa(taxa_names(physeq) %in% p_tax, physeq)
    # use this reduced table to get a frequency table to calculate 
    if (mode == "richness"){
        output <- specnumber(otu_table(red_physeq) %>% as.matrix(), MARGIN = 2)/ntaxa(physeq)
    } else if (mode == "evenness"){
        #output <- diversity(otu_table(red_physeq) %>% as.matrix(), MARGIN = 2, index = "simpson")
        output <- colSums(otu_table(red_physeq))/colSums(otu_table(physeq))
    }
    return(output) 
}

#' Evaluating the richness, eveness and number of traits  
eval_func <- function(physeq, q_sets, sites, t_class){
    if (nrow(q_sets@element) == 0){
        richness <- 0
        evenness <- 0
        n_traits <- 0 
    } else {
        richness <- count_tax(physeq, q_sets, "richness")
        evenness <- count_tax(physeq, q_sets, "evenness")
        n_traits <- count_traits(physeq, q_sets)
    }
    tibble(
        n_traits = n_traits,
        richness = richness,
        evenness = evenness, 
        site = sites, 
        class = t_class
    )
}