from snakemake.utils import min_version 
min_version("6.0")

rule all:
    input: 
        expand("output/coverage/{sequencing}_coverage_{category}.rds", 
            sequencing = ["16s", "wgs"], 
            category=["metabolism", "gram_stain", "sporulation", "motility", 
                        "cell_shape", "substrate", "pathways"])


rule generate_trait_scores: 
    input: 
        set_path=lambda wildcards: "output/sets/set_" + wildcards.category + "_" + 
                        config[wildcards.sequencing] + ".rds"
    params:
        category="{category}",
        sequencing="{sequencing}"
    output:
        "output/coverage/{sequencing}_coverage_{category}.rds"
    script:
        "../../R/coverage_analysis.R"