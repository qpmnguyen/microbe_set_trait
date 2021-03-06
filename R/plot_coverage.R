library(tidyverse)
library(ggthemes)
library(ggsci)
library(MetBrewer)
library(tidytext)
library(ggrepel)
theme_nice <- function() {
    theme_minimal(base_family = "National 2") +
        theme(panel.grid.minor = element_blank(),
              #plot.background = element_rect(fill = "white", color = NA),
              plot.title = element_text(face = "bold"),
              axis.title = element_text(face = "bold"),
              strip.text = element_text(face = "bold", size = rel(0.8), hjust = 0),
              #strip.background = element_rect(fill = "grey80", color = NA),
              legend.title = element_text(face = "bold"))
}
theme_set(theme_nice())

plot_coverage <- function(file_path){
    results <- readRDS(file = file_path)

    results <- results %>% 
        mutate(major_site = map_chr(site, ~{str_split(.x, ":")[[1]][1]})) %>% 
        mutate(minor_site = map_chr(site, ~{str_split(.x, ":")[[1]][2]})) %>% 
        mutate(type = str_to_title(str_replace(type, pattern = "_", replacement = " "))) %>% 
        mutate(minor_site = str_to_title(str_replace(minor_site, pattern = "_", replacement = " ")))
    
    plt <- ggplot(results, 
           aes(x = reorder_within(minor_site, tot_match, type), 
               y = tot_match, fill = major_site)) + 
        geom_bar(stat = "identity") + scale_x_reordered() +
        facet_wrap(~type, scales = "free_y") + 
        coord_flip() + 
        scale_fill_manual(values = met.brewer("Stevens")) + 
        labs(y = "Proportion of identified taxa has an assigned trait", fill = "Major Site", 
             x = "Sites")
    return(plt)
}

plot_16s <- plot_coverage(file.path("output", "coverage", "hmp_16s_cov.rds"))
plot_wgs <- plot_coverage(file.path("output", "coverage", "hmp_wgs_cov.rds"))
plot_wgs

plot_cov_granular <- function(file_path){
    results <- readRDS(file = file_path)
    results <- results %>%
        filter(type %in% c("carbon_substrates", "pathways")) %>%
        mutate(major_site = map_chr(site, ~{str_split(.x, ":")[[1]][1]})) %>% 
        mutate(minor_site = map_chr(site, ~{str_split(.x, ":")[[1]][2]})) %>% 
        mutate(type = str_to_title(str_replace(type, pattern = "_", replacement = " "))) %>% 
        mutate(minor_site = str_to_title(str_replace(minor_site, pattern = "_", replacement = " ")))
    
    results <- results %>% unnest(by_set) %>% mutate(label = "") %>% group_by(type, site) %>% 
        mutate(r = rank(-prop)) %>% ungroup() %>%
        mutate(label = if_else(r <= 3, set_names, ""))
    
    plt <- ggplot(results, aes(y = prop, x = minor_site, col = major_site, label = label)) + 
        geom_text_repel(max.overlaps = 20) +
        geom_jitter() + 
        scale_color_manual(values = met.brewer("Stevens")) +
        facet_wrap(~type, nrow = 2) + coord_flip() + 
        labs(y = "Proportion of identified taxa has an assigned trait", col = "Major Site", 
             x = "Sites")
    return(plt)
}

plot_gr_16s <- plot_cov_granular(file.path("output", "coverage", "hmp_16s_cov.rds"))
plot_gr_wgs <- plot_cov_granular(file.path("output", "coverage", "hmp_wgs_cov.rds"))
plot_gr_wgs
    
