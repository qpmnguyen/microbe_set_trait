from snakemake.utils import min_version
min_version("6.0")
#DATASETS=["ibd_16s", "crc_16s"]

rule all: 
    input:
        expand("output/sequence_process_16s/{dataset}/exports/taxonomy.biom", dataset=["ibd_16s", "crc_16s"]),
        expand("output/sequence_process_16s/{dataset}/exports/feature_table.biom", dataset=["ibd_16s", "crc_16s"]),
        expand("output/sequence_process_16s/{dataset}/exports/dna_sequences.fasta", dataset=["ibd_16s", "crc_16s"])



rule import_manifest:
    input: 
        lambda wildcards: config["data"][wildcards.dataset]["manifest"]
    conda: 
        "../env/qiime2-2022.2-py38-linux-conda.yml"
    output:
        temp("output/sequence_process_16s/{dataset}/sequence-demux.qza")
    params:
        input_type=lambda wildcards: config["data"][wildcards.dataset]["type"],
        input_format=lambda wildcards: config["data"][wildcards.dataset]["format"]
    shell:
        """
        qiime tools import \
            --type {params.input_type} \
            --input-path {input} \
            --output-path {output} \
            --input-format {params.input_format}
        """

rule denoise: 
    input:
        "output/sequence_process_16s/{dataset}/sequence-demux.qza"
    conda: 
        "../env/qiime2-2022.2-py38-linux-conda.yml"
    output:
        table="output/sequence_process_16s/{dataset}/feature_table.qza",
        sequences="output/sequence_process_16s/{dataset}/dna_sequences.qza", 
        stats="output/sequence_process_16s/{dataset}/denoise_stats.qza"
    threads: 10
    params: 
        quality=lambda wildcards: config["data"][wildcards.dataset]["params"],
        cmd=lambda wildcards: config["data"][wildcards.dataset]["cmd"]
    shell:
        """
        qiime dada2 {params.cmd} \
            --i-demultiplexed-seqs {input} \
            {params.quality} \
            --p-trunc-q 2 \
            --p-pooling-method 'pseudo' \
            --p-chimera-method 'pooled' \
            --p-n-threads {threads} \
            --o-table {output.table} \
            --o-representative-sequences {output.sequences} \
            --o-denoising-stats {output.stats}
        """

rule download_classifier: 
    output:
        "large_files/silva-138-99-nb-weighted-classifier.qza"
    shell:
        """
        wget https://data.qiime2.org/2022.2/common/silva-138-99-nb-weighted-classifier.qza \
            -O {output}
        """
        
rule taxonomic_classification:
    input:
        classifier="large_files/silva-138-99-nb-weighted-classifier.qza",
        reads="output/sequence_process_16s/{dataset}/dna_sequences.qza"
    conda: 
        "../env/qiime2-2022.2-py38-linux-conda.yml"
    output:
        "output/sequence_process_16s/{dataset}/taxonomy.qza"
    shell:
        """
        qiime feature-classifier classify-sklearn \
            --i-classifier {input.classifier} \
            --i-reads {input.reads} \
            --o-classification {output}
        """
    
rule export:
    input:
        feature_table="output/sequence_process_16s/{dataset}/feature_table.qza",
        taxonomy="output/sequence_process_16s/{dataset}/taxonomy.qza",
        sequences="output/sequence_process_16s/{dataset}/dna_sequences.qza"
    conda: 
        "../env/qiime2-2022.2-py38-linux-conda.yml"
    output:
        outtax="output/sequence_process_16s/{dataset}/exports/taxonomy.tsv", 
        outtable="output/sequence_process_16s/{dataset}/exports/feature-table.biom",
        outseq="output/sequence_process_16s/{dataset}/exports/dna-sequences.fasta"
    shell:
        """
        qiime tools export \
            --input-path {input.feature_table} \
            --output-path output/sequence_process_16s/{wildcards.dataset}/exports       
        qiime tools export \
            --input-path {input.taxonomy} \
            --output-path output/sequence_process_16s/{wildcards.dataset}/exports
        qiime tools export \
            --input-path {input.sequences} \
            --output-path output/sequence_process_16s/{wildcards.dataset}/exports
        """
        
