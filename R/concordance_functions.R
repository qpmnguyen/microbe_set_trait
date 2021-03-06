library(phyloseq)
library(tidyverse)
library(glue)

t_test <- function(trait_data, met, outcome_var, pos_class, ret = "p.value", correction = "BH"){
    outcome <- met %>% pull(outcome_var)
    x_dat <- trait_data[which(outcome == pos_class),]
    y_dat <- trait_data[-which(outcome == pos_class),]
    results <- map_dbl(seq_len(ncol(x_dat)), ~{
        mod <- t.test(x_dat[,.x], y_dat[,.y])
        mod[[ret]]
    })
    names(results) <- colnames(x_dat)
    if (ret == "p.value"){
        results <- p.adjust(results, method = correction)
    }
    return(results)
}


# pathways first 
#' @param trait_data The trait csv data set as data.frame
#' @param path_data The path data set as a data.frame
#' @param met The metadata annotation as data.frame
#' @param ref_df The set-to-reference data frame where set is the set label and db is the name in MetaCyc
#' @param o_var Outcome variable
#' @param p_class Positive class 
#' @param metadata The list of pathways by annotation
#' @param annotation Any additional annotation to attach to the last column
assess_traits <- function(trait_data, path_data, met, ref_df, o_var, p_class, metadata, annotation){
    ref_ids <- ref_df %>% pull(set)
    trait_data <- trait_data %>% dplyr::select(which(colnames(trait_data) %in% ref_ids)) 
    # get significant traits 
    sig_traits <- t_test(trait_data, met, outcome_var = o_var, pos_class = p_class)
    sig_traits <- names(sig_traits)[sig_traits <= 0.05]
    # sig_db 
    sig_db <- ref_df %>% filter(set %in% sig_traits) %>% pull(db)

    
    # extract pathways 
    subset_metadata <- metadata[sig_db]
    retrieve_pathways <- map(subset_metadata, ~which(colnames(path_data) %in% .x))
    
    results <- imap(retrieve_pathways, ~{
        if (length(.x) == 0){
            out <- tibble(trait = .y, prop = NA_real_, size = NA_real_)
        } else {
            # per model 
            p_values <- vector(length = length(.x))
            t_out <- ifelse(met %>% pull(o_var) == p_class,1,0)
            for (i in seq_along(.x)){
                t_vec <- path_data[,.x[i]]
                t_vec <- asin(sqrt(abs(t_vec/sum(t_vec))))
                mod <- lm(t_vec ~ as.factor(t_out))
                p_values[i] <- coef(summary(mod))[,4][-1]
            }
            p_values <- p.adjust(p_values, method = "BH")
            prop <- length(which(p_values <= 0.05))/length(p_values)
            out <- tibble(trait = .y, prop = prop, size = length(.x))
        } 
        return(out)
    })
    results <- do.call(bind_rows, results)
    results <- results %>% mutate(type = rep(annotation, nrow(results)))
    return(results)
}

assess_correlation <- function(trait_data, path_data, ref_df, metadata){
    ref_ids <- ref_df %>% pull(set) %>% unique()
    
    common <- intersect(colnames(trait_data), ref_ids)

    results <- imap(common, ~{
        t_vec <- trait_data %>% pull(.x)
        db_intersect <- ref_df %>% filter(set == .x) %>% pull(db)
        db_intersect <- metadata[[db_intersect]]
        p_df <- path_data %>% select(any_of(db_intersect))
        if (ncol(p_df) == 0){
            out <- tibble(trait = .x, corr = NA_real_, size = NA_real_, pathway = NA_character_)
        } else {
            corr <- vector(length = ncol(p_df))
            pval <- vector(length = ncol(p_df))
            for (i in seq_len(ncol(p_df))){
                mod <- cor.test(x = t_vec, y = p_df[,i], method = "spearman")
                corr[i] <- mod$estimate
                pval[i] <- mod$p.value
            }
            pval <- p.adjust(pval, method = "BH")
            out <- tibble(trait = .x, corr = corr, size = ncol(p_df), pathway = colnames(p_df), pval = pval)
        }
        return(out)
    })
    results <- do.call(bind_rows, results)
    results <- results %>% left_join(annotation_df) %>% drop_na(corr)
    results <- results %>% mutate(pnames = str_wrap(paste(pathway, annotation, sep = ": "), width = 20))
    
    corr_mat <- results %>% select(trait, corr, pnames) %>% 
        pivot_wider(names_from = "pnames", values_from = "corr") %>% 
        column_to_rownames("trait") %>% as.matrix()
    
    p_mat <- results %>% select(trait, pval, pnames) %>% 
        pivot_wider(names_from = "pnames", values_from = "pval") %>% 
        column_to_rownames("trait") %>% as.matrix()
    return(
        list(
            results = results, 
            corr_mat = corr_mat, 
            p_mat = p_mat
        )
    )
}



# old Maaslin2 code 

# else {
#     subset_path <- path_data %>% dplyr::select(.x)
#     invisible(capture.output(mod <- Maaslin2(
#         input_data = subset_path,
#         input_metadata = met,
#         output = tempdir(),
#         normalization = "TSS",
#         transform = "AST",
#         analysis_method = "LM", fixed_effects = "study_condition",
#         reference = "study_condition, control",
#         correction = "BH")))
#     pvals <- mod$results %>% pull(qval)
#     prop <- which(pvals <= 0.05) %>% length() %>% `/`(.,nrow(mod$results))
#     out <- tibble(trait = .y, prop = prop, size = length(.x))
# }

