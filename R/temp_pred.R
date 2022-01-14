source("R/functions_pred.R")


wgs <- get_ihmp(type = "traits")
amplicon <- get_gevers()

wgs_gsva <- fit_gsva(wgs, t_rank = "species")
amplicon_gsva <- fit_gsva(amplicon, t_rank = "genus")


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
    best_rf <- rand_forest(min_n = !!param, trees = 3000) %>% 
        set_engine("ranger", num_threads = !!n_threads) %>%
        set_mode("classification")
    
    new_wkflow <- wkflow %>% update_model(best_rf)
    return(new_wkflow %>% last_fit(split))
}

wgs_fit <- fit_final(wgs_gsva)


