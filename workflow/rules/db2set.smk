from snakemake.utils import min_version 
min_version("6.0")

rule all: 
    input:
        expand("output/sets/set_{category}_{level}.rds", 
        category=["metabolism", "gram_stain", "sporulation", "motility", "cell_shape", "substrate", "pathways"], 
        level=["species", "genus"])

rule retrieve_sets:
    input: 
        database="output/databases/madin_proc.rds",
    params:
        category="{category}",
        level="{level}",
        filt_thresh=0.95
    output:
        "output/sets/set_{category}_{level}.rds"
    script:
        "../../R/db_2_set.R"
