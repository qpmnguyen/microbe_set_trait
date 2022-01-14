library(tidymodels)
library(GSVA)
library(glue)
library(curatedMetagenomicData)
library(BiocSet)
library(themis)
library(compositions)
requireNamespace("speedyseq")
requireNamespace("randomForest")
requireNamespace("mia")
source("R/utils.R")
source("R/gevers_preprocess.R")
source("R/functions_coverage.R")

# DATA LOADING #### 
get_gevers <- function(){
    physeq <- load_gevers() 
    physeq <- attach_ncbi_std(physeq, t_rank = "Genus")
    physeq <- phyloseq::transform_sample_counts(physeq, function(x) x + 1)
    return(physeq)
}

get_picrust2 <- function(metadata, type){
    type <- match.arg(type, c("pathway", "enzyme"))
    if (type == "pathway"){
        tbl <- fread(file.path("picrust2_output", "path_abun_unstrat_desc.tsv.gz"))
        pred_var <- "pathway"
    } else {
        tbl <- fread(file.path("picrust2_output", "pred_metagenome_unstrat_desc.tsv.gz"))
        pred_var <- "function"
    } 
    df <- tbl %>% as.data.frame() %>% select(-description) %>% 
        column_to_rownames(var = pred_var) %>% t() 
    
    df <- clr(acomp(df + 1)) %>% as.data.frame() %>% 
        rownames_to_column("sample_id") %>% 
        left_join(metadata)
    
    return(df)
}

get_ihmp <- function(type) {
    type <- match.arg(type, c("pathway", "function", "traits"))
    if (type == "pathway"){
        physeq <- curatedMetagenomicData(pattern = "HMP_2019_ibdmdb.pathway_abundance", 
                               dryrun = FALSE, rownames = "NCBI")[[1]]
        # remove stratified results 
        df <- assay(physeq) %>% as.data.frame() %>% rownames_to_column(var = "pathway") %>% 
            filter(str_detect(pathway, "\\|", negate = TRUE)) %>% 
            column_to_rownames(var = "pathway") %>% t() 
        # perform centered log ratio transformation
        df <- clr(acomp(df + 10e-5)) %>% as.data.frame() %>% 
            rownames_to_column("sample_id") %>% as_tibble()
        
        metadata <- colData(physeq) %>% as.data.frame() %>% 
            rownames_to_column(var = "sample_id") %>% 
            mutate(diagnosis = disease) %>% select(sample_id, diagnosis)  
        df <- left_join(df, metadata)
    } else if (type == "function"){
        physeq <- curatedMetagenomicData(pattern = "HMP_2019_ibdmdb.gene_families", dryrun = FALSE)
    } else if (type == "traits"){
        physeq <- curatedMetagenomicData(pattern = "HMP_2019_ibdmdb.relative_abundance", 
                               dryrun = FALSE, rownames = "NCBI")[[1]] %>% 
            mia::makePhyloseqFromTreeSummarizedExperiment(abund_values = "relative_abundance")
        physeq <- attach_ncbi_metaphlan(physeq)
    }
    # return output
    if (type == "traits"){
        return(physeq)
    } else {
        return(df)
    }
}

get_metadata <- function(physeq, class){
    class <- match.arg(class, c("16s", "wgs"))
    if (class == "16s"){
        s_data <- sample_data(physeq) %>% as("data.frame") %>% 
            rownames_to_column(var = "sample_id") %>% select(sample_id, diagnosis)
    } else {
        s_data <- sample_data(physeq) %>% as("data.frame") %>% rownames_to_column(var = "sample_id") %>% 
            mutate(diagnosis = disease) %>% 
            select(sample_id, diagnosis)
    } 
    return(s_data)
}

# PRODUCING SETS ####
# this returns a list as a set with list of elements 
match_sets <- function(physeq, sets, out_type = "list"){
    physeq <- subset_taxa(physeq, !is.na(ncbiids))
    physeq <- speedyseq::tax_glom(physeq, taxrank = "ncbiids")
    taxtab <- as(tax_table(physeq),"data.frame")
    list_set <- as(sets, "list")
    list_set <- map(list_set, ~{
        rownames(taxtab)[which(taxtab$ncbiids %in% .x)]
    })
    list_set <- Filter(length, list_set)
    if (out_type == "list"){
        return(list_set)
    } else {
        return(BiocSet::BiocSet(list_set))
    }
}

# this function creates a union of all the sets from all types
unify_sets <- function(physeq, t_rank){
    t_rank <- match.arg(t_rank, c("genus", "species"))
    type <- c("carbon_substrates", "cell_shape", "gram_stain", "motility", 
             "pathways", "range_salinity", "range_tmp", "sporulation")
    set_list <- map(type, ~{
        sets <- readRDS(file = file.path("output", "sets", glue("madin_{type}_{t_rank}.rds", 
                                                        type = .x, t_rank = t_rank)))
        tmp <- as(sets, "list")
        names(tmp) <- paste(.x, names(tmp), sep = ":")
        return(BiocSet::BiocSet(tmp))
    })
    cov <- imap(set_list, ~calculate_coverage(physeq, sets = .x, type = type[.y], site = NA_character_))
    cov <- do.call(bind_rows, cov)
    
    big_set <- Reduce(generics::union, set_list)
    output <- list(
        cov = cov, 
        set = big_set
    )
    return(output)
}



# FITTING MODELS #### 
# fit gsva for the model 
fit_gsva <- function(physeq, t_rank){
    big_set <- unify_sets(physeq, t_rank = t_rank)$set
    filt_set <- match_sets(physeq, big_set)
    otu <- as(otu_table(physeq), "matrix")
    if (t_rank == "genus"){
        otu <- t(otu)
        kernel <- "Poisson"
        class <- "16s"
    } else {
        kernel <- "Gaussian"
        class <- "wgs"
    }
    enrich_mat <- t(gsva(expr = otu, gset.idx.list = filt_set, 
                         method = "gsva", kcdf = "Gaussian", min.sz = 2)) %>% 
        as.data.frame() %>% rownames_to_column(var = "sample_id")
    metadata <- get_metadata(physeq, class = class) 
    enrich_mat <- left_join(enrich_mat, metadata)
    return(enrich_mat)
}

define_modelflow <- function(data, n_threads){
    proc_recipe <- recipe(diagnosis ~ ., data = data) %>% 
        update_role(sample_id, new_role = "sample ID") %>%
        step_smote(diagnosis) %>% 
        step_normalize(all_numeric_predictors())
    
    rf_mod <- rand_forest(min_n = tune(), trees = 3000) %>% 
        set_engine("ranger", num.threads = !!n_threads) %>%
        set_mode("classification")
    
    wkflow <- workflow() %>%
        add_model(rf_mod) %>%
        add_recipe(proc_recipe)
    return(wkflow)
}

complete_fit <- function(df, grid_size, init_wkflow, annotate = NULL, n_threads){
    nested_eval <- nested_cv(df, outside = vfold_cv(v = 10), inside = vfold_cv(v = 5))
    print("Starting fit...")
    outer_fit <- pmap_dfr(nested_eval, function(splits, id, inner_resamples){
        print("Starting tune...")
        tune <- fit_models(resamp = inner_resamples, 
                           wkflow = init_wkflow, grid_size = grid_size, 
                           n_threads = n_threads)
        print("Ending tune ...")
        tune$adj_wkflow %>% last_fit(splits) %>% collect_metrics() %>% 
            filter(.metric == "roc_auc")
    })
    print("Finished fitting...")
    if (!is.null(annotate)){
        outer_fit <- outer_fit %>% mutate(annotate = annotate)
    }
    return(outer_fit)
}

# fit a model based on training data
fit_models <- function(resamp, wkflow, grid_size, n_threads){
    #val_splits <- vfold_cv(train_dat, v = 5)
    #TODO: Review how to evaluate binary classification models
    results <- wkflow %>% tune_grid(resamples = resamp, grid = grid_size, 
                                    control = control_grid(save_pred = TRUE),
                                    metrics = metric_set(roc_auc))
    best_mod <- results %>% select_best(metric = "roc_auc")
    internal_pred <- results %>% collect_predictions(parameters = best_mod) 
    
    param <- best_mod$min_n
    best_mod <- rand_forest(min_n = !!param, trees = 3000) %>% 
        set_engine("ranger", num.threads = !!n_threads) %>%
        set_mode("classification")
    return(list(
        overall = results %>% collect_metrics(),
        model = best_mod,
        internal = internal_pred,
        adj_wkflow = wkflow %>% update_model(best_mod)
    ))
}

