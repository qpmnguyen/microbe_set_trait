from snakemake.utils import min_version
min_version("6.0")

rule all: 


rule import_manifest:
    input: 
        "metadata/crc_qiime2_metadata.tsv"
    conda: 
        "../env/qiime2-2022.2-py38-linux-conda.yml"
    output:
        "output/sequence_process_16s/crc_16s/paired-end-dada2.qza"
    shell:
        """
        qiime tools import \\
            --type 'SampleData[PairedEndSequencesWithQuality]' \\
            --input-path {input}
            --output-path {output}
            --input-format PairedEndFastqManifestPhred33V2
        """

rule dada2_run: 
    conda: 
        "../env/qiime2-2022.2-py38-linux-conda.yml"

rule taxonomic_classification:
    conda: 
        "../env/qiime2-2022.2-py38-linux-conda.yml"
    
rule export:
    conda: 
        "../env/qiime2-2022.2-py38-linux-conda.yml"