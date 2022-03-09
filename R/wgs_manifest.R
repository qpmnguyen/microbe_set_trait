# Downloading wgs data sets for prediction  
# Quang Nguyen 
# Last revision: 2022-03-09
library(curatedMetagenomicData)
library(tidyverse)

metadata <- as_tibble(sampleMetadata)
s_names <- metadata %>% filter(study_condition %in% c("CRC", "T2D", "IBD")) %>% pull(study_name) %>% unique()
samples <- metadata %>% filter(study_name %in% s_names, study_condition %in% c("CRC", "T2D", "IBD", "control"))


control_only <- metadata %>% filter(study_condition == "control")
