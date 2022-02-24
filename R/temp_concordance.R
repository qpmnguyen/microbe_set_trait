library(tidyverse)
library(fgsea)
library(ANCOMBC)
library(gt)
wgs <- get_ihmp(type = "traits")
wgs_path <- get_ihmp(type = "pathway")


wgs_path

big_set <- unify_sets(wgs, t_rank = "species")$set
filt_set <- match_sets(wgs, big_set)

mod_wgs <- ancombc(wgs, formula = "disease", p_adj_method = "holm", 
        group = "disease", neg_lb = TRUE) 

coef <- mod_wgs$res$beta

input <- coef$diseaseIBD
names(input) <- rownames(coef)
input <- input[order(input, decreasing = FALSE)]
wgs_enrich <- fgsea(pathways = filt_set, stats = input, minSize = 2)

top_up <- wgs_enrich[ES > 0][head(order(pval), n = 10),pathway]
top_down <- wgs_enrich[ES < 0][head(order(pval), n = 10), pathway]

plotGseaTable(filt_set[c(top_up, rev(top_down))], input, wgs_enrich, gseaParam = 0.5)
gt

wgs_enrich[pathway %in% c(top_up, rev(top_down))] %>% 
    mutate(updown = if_else(pathway %in% top_up, "Enriched", "Depleted")) %>%
    group_by(updown) %>%
    arrange(ES, pval, .by_group = TRUE) %>% 
    select(pathway, pval, padj, NES, updown) %>%
    gt(rowname_col = "pathway", groupname_col = "updown") %>%
    tab_header(title = md("Enrichment analysis using *ANCOM-BC* and *fgsea*"), 
               subtitle = "Top 10 most enriched/depleted traits") %>%
    fmt_number(columns = c(pval, padj, NES)) %>%
    cols_label(pval = "Raw p-values", padj = "Adjusted p-values", NES = "Effect Size") %>%
    tab_options(table.font.names = "National 2") %>%
    tab_style(
        style = list(cell_fill(color = "salmon", alpha = 0.5), cell_text(weight = "bold")),
        locations = cells_row_groups()
    ) %>% 
    tab_style(
        style = list(cell_fill(color = "grey90", alpha = 0.5)),
        locations = cells_body(columns = "padj")
    ) %>% gtsave(filename = "temp_wgs_trait.png")
    

