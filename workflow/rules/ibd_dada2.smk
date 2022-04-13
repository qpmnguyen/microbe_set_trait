from snakemake.utils import min_version
min_version("6.0")

rule all: 


rule import_manifest:
    input: 
        "metadata/crc_16s/
    conda: 
        "../env/qiime2-2022.2-py38-linux-conda.yml"
    output:
    shell:

rule dada2_run: 
    conda: 
        "../env/qiime2-2022.2-py38-linux-conda.yml"

rule taxonomic_classification:
    conda: 
        "../env/qiime2-2022.2-py38-linux-conda.yml"
    
rule export:
    conda: 
        "../env/qiime2-2022.2-py38-linux-conda.yml"