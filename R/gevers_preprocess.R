library(tidyverse)


mac_path <- "~/research/data/ena_gevers_et_al/master/"

files <- list.files(mac_path)

sample_names <- map_chr(str_split(files, pattern = ".fastq.gz"), ~{ .x[1]})
gevers_metadata <- read.table(file = "metadata/gevers_metadata.txt", sep = "\t", header = TRUE) 
sites <- unique(gevers_metadata$biopsy_location)

# create a breakdown of diagnosis per site
breakdown <- map(sites, ~{
  gevers_metadata %>% filter(biopsy_location == .x) %>% 
    group_by(diagnosis) %>% count()
})
names(breakdown) <- sites

# terminal illeum has the best amount of data available, with balanced evaluation between controls (no IBD) 
# and crohn's disease. 
sample_metadata <- gevers_metadata %>% filter(biopsy_location == "Terminal ileum") %>% 
  select(sample_name, diagnosis) %>% 
  filter(diagnosis %in% c("no", "CD")) 
common <- intersect(sample_names, sample_metadata$sample_name)

# create a new directory for all matching samples 
dir.create("~/research/data/ena_gevers_et_al/input_data")
for (i in list.files(mac_path, full.names = TRUE)){
  names <- strsplit(strsplit(i,"//")[[1]][2], ".fastq.gz", fixed = TRUE)[[1]]
  if (names %in% common){
    file.copy(i, "~/research/data/ena_gevers_et_al/input_data")
  }
}

