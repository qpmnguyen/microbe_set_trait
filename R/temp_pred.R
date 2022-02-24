source("R/functions_pred.R")
library(tidyverse)
library(ggthemes)
library(ggsci)
library(MetBrewer)
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

wgs <- get_ihmp(type = "traits")
amplicon <- get_gevers()

wgs_gsva <- fit_gsva(wgs, t_rank = "species")
amplicon_gsva <- fit_gsva(amplicon, t_rank = "genus")

wgs_path <- get_ihmp(type = "pathway")
amplicon_path <- get_picrust2(metadata = get_metadata(amplicon, class = "16s"), 
                              type = "pathway")



wkflow <- define_modelflow(data = wgs_gsva, n_threads = 2)

fit_final <- function(data){
    wkflow <- define_modelflow(data = data, n_threads = 2)
    split <- initial_split(data, prop = 0.8)
    train_dat <- training(split)
    test_dat <- testing(split)
    
    val_split <- vfold_cv(train_dat, v = 5)
    
    results <- wkflow %>% tune_grid(resamples = val_split, grid = 25, 
              control = control_grid(save_pred = TRUE),
              metrics = metric_set(roc_auc))
    
    best_mod <- results %>% select_best(metric = "roc_auc")
    param <- best_mod$min_n
    best_rf <- rand_forest(min_n = !!param, trees = 2000) %>% 
        set_engine("ranger", num.threads = 2, importance = "impurity") %>%
        set_mode("classification")
    
    new_wkflow <- wkflow %>% update_model(best_rf)
    return(new_wkflow %>% last_fit(split))
}

wgs_fit <- fit_final(wgs_gsva)
wgs_vip_features <- wgs_fit %>% pluck(".workflow",1) %>% extract_fit_parsnip() %>% vip(num_features = 15)
wgs_vip_features + geom_bar(fill = "steelblue", stat = "identity") + 
    theme(axis.text.y = element_text(size = 10)) + 
    labs(title = "Top 15 most important trait features for iHMP data set", 
         subtitle = paste("AUROC:", round(wgs_fit %>% collect_metrics() %>% 
                                              filter(.metric == "roc_auc") %>% 
                                              pull(.estimate),3)))

amplicon_fit <- fit_final(amplicon_gsva)
amplicon_vip_features <- amplicon_fit %>% pluck(".workflow", 1) %>% extract_fit_parsnip() %>% vip(num_features = 15)
amplicon_vip_features + geom_bar(fill = "salmon", stat = "identity") + 
    theme(axis.text.y = element_text(size = 10)) + 
    labs(title = "Top 15 most important features for the Gevers et al. data set", 
         subtitle = paste("AUROC:", round(amplicon_fit %>% collect_metrics() %>% 
                                              filter(.metric == "roc_auc") %>% 
                                              pull(.estimate),3)))

wgs_path_fit <- fit_final(wgs_path)
amplicon_path_fit <- fit_final(amplicon_path)

wgs_path_vip <- wgs_path_fit %>% pluck(".workflow", 1) %>% extract_fit_parsnip() %>% vip(num_features = 15)
wgs_path_vip + geom_bar(fill = "salmon", stat = "identity") + 
    theme(axis.text.y = element_text(size = 10)) + 
    labs(title = "Top 15 most important pathway features for the iHMP data set", 
         subtitle = paste("AUROC:", round(wgs_path_fit %>% collect_metrics() %>% 
                                              filter(.metric == "roc_auc") %>% 
                                              pull(.estimate),3)))

df_result <- bind_rows(
    amplicon_path_fit %>% collect_metrics() %>% 
        mutate(type = "16s_path"),
    wgs_path_fit %>% collect_metrics() %>% 
        mutate(type = "wgs_path"),
    wgs_fit %>% collect_metrics() %>% 
        mutate(type = "wgs_trait"),
    amplicon_fit %>% collect_metrics() %>%
        mutate(type = "16s_trait")
)

custom_colors <- c("#FA8072", "#4682b4")
names(custom_colors) <- c("16S", "WGS")

ggplot(df_result %>% filter(.metric == "roc_auc") %>% 
           mutate(data_type = map_chr(type, ~str_split(.x, "_")[[1]][1])) %>%
           mutate(data_type = str_to_upper(data_type)) %>% 
           mutate(source = map_chr(type, ~str_split(.x, "_")[[1]][2])) %>%
           mutate(source = if_else(source == "path", "Pathway Abundances", "Trait Enrichment Scores")), 
       aes(y = .estimate, x = source, fill = source)) + 
    scale_fill_manual(values = c('salmon', "steelblue")) +
    geom_bar(stat = "identity") + 
    facet_wrap(~data_type) +
    labs(y = "Test set AUROC", x = "Data Source") +
    theme(legend.position = "None", 
          strip.text.x = element_text(size = 13))



