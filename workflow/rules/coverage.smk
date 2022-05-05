from snakemake.utils import min_version 
min_version("6.0")

rule all:
    input: 
        expand("output/coverage/{dset}_coverage_{category}.rds", 
            dset = ["hmp_16s", "hmp_wgs", "ibd_16s", "crc_16s", "ibd_wgs", "ibd_wgs"], 
            category=["metabolism", "gram_stain", "sporulation", "motility", 
                        "cell_shape", "substrate", "pathways"])


rule generate_trait_scores: 
    input: 
        set_path=lambda wildcards: "output/sets/set_" + wildcards.category + "_" + 
                        config[wildcards.dset] + ".rds"
    params:
        category="{category}",
        dset="{dset}"
    output:
        "output/coverage/{dset}_coverage_{category}.rds"
    script:
        "../../R/coverage_analysis.R"