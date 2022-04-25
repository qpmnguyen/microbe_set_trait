from snakemake.utils import min_version 
min_version("6.0")
DATASETS = ["ibd_16s", "crc_16s"]


rule all: 
    input:
        expand("output/sets/set_{class}_{level}.rds", datasets = DATASETS)
