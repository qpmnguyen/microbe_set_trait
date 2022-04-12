from snakemake.utils import min_version
min_version("6.0")

rule all:
    input: 
        "notebooks/tester.ipynb"
    output:
        "notebooks/docs/tester.md"
    shell:
        "jupyter nbconvert --execute {input} --ExecutePreprocessor.kernel_name=ir --output docs/tester.md --to markdown"

    

