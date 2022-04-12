from snakemake.utils import min_version
min_version("6.0")

rule all:
    conda:
        "../env/db_prep.yml"
    output: 
        "../../output/databases/db_merged.csv"
    log:
        notebook="../../logs/notebooks/db_prep.r.ipynb"
    notebook:
        "../../notebooks/db_prep.r.ipynb"

