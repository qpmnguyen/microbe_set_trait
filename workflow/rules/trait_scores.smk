from snakemake.utils import min_version 
min_version("6.0")

rule all:
    input: 
        expand("data/trait_{condition}_{sequencing}_feat.csv", condition = ["crc", "ibd"], 
            sequencing = ["16s", "wgs"])


rule generate_trait_scores: 
    input: 
        "data/pred_relabun_{condition}_{sequencing}.rds"
    params:
        condition="{condition}",
        sequencing="{sequencing}"
    output:
        "data/trait_{condition}_{sequencing}_feat.csv"
    script:
        "../../R/generate_scores.R"