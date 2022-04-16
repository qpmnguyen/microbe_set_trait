from snakemake.utils import min_version 
min_version("6.0")
DATASETS = ["ibd_16s", "crc_16s"]

rule all:
    input:
        expand("output/picrust2/{datasets}/pathways_out/path_abun_unstrat_desc.tsv.gz", datasets = DATASETS)


rule place_seqs:
    input: 
        "output/sequence_process_16s/{datasets}/exports/dna-sequences.fasta"
    conda: 
        "../env/picrust_env.yml"
    output:
        "output/picrust2/{datasets}/placed_seqs.tre"
    threads: 4
    shell:
        "place_seqs.py -s {input} -o {output} -p {threads} --verbose -t epa-ng"

rule predict_ec:
    input:
        "output/picrust2/{datasets}/placed_seqs.tre"
    conda: 
        "../env/picrust_env.yml"
    output:
        "output/picrust2/{datasets}/EC_predicted.tsv.gz"
    threads: 4
    shell:
        "hsp.py -i EC -t {input} -o {output} -p {threads}"

rule predict_copy:
    input: 
        "output/picrust2/{datasets}/placed_seqs.tre"
    conda: 
        "../env/picrust_env.yml"
    output:
        "output/picrust2/{datasets}/marker_predicted_and_nsti.tsv.gz"
    threads: 4
    shell:
        "hsp.py -i 16S -t {input} -o {output} -p {threads}"

rule metagenome_prediction:
    input:
        biom="output/sequence_process_16s/{datasets}/exports/feature-table.biom",
        marker="output/picrust2/{datasets}/marker_predicted_and_nsti.tsv.gz",
        ec="output/picrust2/{datasets}/EC_predicted.tsv.gz"
    conda: 
        "../env/picrust_env.yml"
    output:
        "output/picrust2/{datasets}/EC_metagenome_out/pred_metagenome_unstrat.tsv.gz",
        "output/picrust2/{datasets}/EC_metagenome_out/seqtab_norm.tsv.gz",
        "output/picrust2/{datasets}/EC_metagenome_out/weighted_nsti.tsv.gz"
    shell:
        """
        metagenome_pipeline.py -i {input.biom} -m {input.marker} \\
            -f {input.ec} \\
            -o output/picrust2/{wildcards.datasets}/EC_metagenome_out
        """

rule predict_pathway:
    input:
        "output/picrust2/{datasets}/EC_metagenome_out/pred_metagenome_unstrat.tsv.gz"
    conda: 
        "../env/picrust_env.yml"
    output:
        "output/picrust2/{datasets}/pathways_out/path_abun_unstrat.tsv.gz"
    threads: 4
    shell:
        "pathway_pipeline.py -i {input} -p output/picrust2/{wildcards.datasets}/pathways_out -p {threads}"

rule add_desc:
    input:
        "output/picrust2/{datasets}/pathways_out/path_abun_unstrat.tsv.gz"
    conda: 
        "../env/picrust_env.yml"
    output:
        "output/picrust2/{datasets}/pathways_out/path_abun_unstrat_desc.tsv.gz"
    shell:
        "add_descriptions.py -i {input} -o {output} -m METACYC"
