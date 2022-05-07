from snakemake.utils import min_version 
min_version("6.0")

CONDITION=["ibd", "crc"]
INPUT_TYPE=["picrust2", "pathways"]
SEQUENCING=["16s", "wgs"]


rule all: 
    input:
        expand("output/pred/fimportance_{condition}_{sequencing}_{input_type}.csv", zip,
        condition=["ibd", "ibd", "crc", "crc", "ibd", "ibd", "crc", "crc"], 
        sequencing=["16s", "wgs", "16s", "wgs", "16s", "wgs", "16s", "wgs"],
        input_type=["picrust2", "pathways", "picrust2", "pathways", "trait", "trait", "trait", "trait"])

rule feature_importance:
    conda:
        "../env/predictive_analysis.yml"
    input: 
        feat_path=lambda wildcards: config[wildcards.condition + "_" +
                                           wildcards.sequencing + "_" + 
                                           wildcards.input_type]["feat_path"],
        meta_path=lambda wildcards: config[wildcards.condition + "_" +
                                           wildcards.sequencing + "_" +
                                           wildcards.input_type]["meta_path"]
    params:
        outcome_label=lambda wildcards: config[wildcards.condition + "_" +
                                               wildcards.sequencing + "_" +
                                               wildcards.input_type]["outcome_label"],
        pos_class=lambda wildcards: config[wildcards.condition + "_" +
                                               wildcards.sequencing + "_" +
                                               wildcards.input_type]["pos_class"],
        condition="{condition}",
        sequencing="{sequencing}",
        input_type="{input_type}",
        clr_trans=lambda wildcards: config[wildcards.condition + "_" +
                                               wildcards.sequencing + "_" +
                                               wildcards.input_type]["clr_trans"]
    threads: 1
    output:
        "output/pred/fimportance_{condition}_{sequencing}_{input_type}.csv"
    script:
        "../../python/feature_importance.py"