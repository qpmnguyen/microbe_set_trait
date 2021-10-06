#!/usr/bin/env nextflow

params.threads = 1
biom_file = Channel.fromPath("data/gevers_biom.biom")
seqs_file = Channel.fromPath("gevers_seq.fa")

process place_seqs {
    publishDir 'picrust2_tmp/', mode: 'copy'
    input:
    file seqs from seqs_file
    output: 
    file "placed_seqs.tre" into placed_tree
    script: 
        """
        place_seqs.py -s $seqs -o placed_seqs.tre -t epa-ng -p ${params.threads}
        """
    stub: 
        """
        touch placed_seqs.tre
        """
}

process hidden_state {
    publishDir 'picrust2_tmp/', mode: 'copy'
    input: 
    file tree from placed_tree
    output:
    file "pred_traits_marker.tsv.gz" into pred_marker 
    file "pred_traits_EC.tsv.gz" into pred_EC 
    script:
        """
        hsp.py -t $tree -i 16S -o pred_traits_marker.tsv.gz -m mp --check -p ${params.threads} --calculate_NSTI
        hsp.py -t $tree -i EC -o pred_traits_EC.tsv.gz -m mp --check -p ${params.threads}
        """
    stub: 
        """
        touch pred_traits_marker.tsv.gz
        touch pred_traits_EC.tsv.gz
        """
}

process metagene_pred {
    publishDir 'picrust2_tmp', mode: 'copy'
    input: 
    file biom from biom_file 
    file marker from pred_marker
    file enzyme from pred_EC 
    output: 
    file 'EC_metagenome_out/pred_metagenome_strat.tsv.gz' into output_EC, output_EC_desc 
    script: 
        """
        metagenome_pipeline.py -i $biom \
        -m ${pred_marker} \
        -f ${pred_EC} \
        -o EC_metagenome_out \
        --strat_out \
        """
    stub:
        """
        mkdir -p EC_metagenome_out
        touch EC_metagenome_out/pred_metagenome_strat.tsv.gz
        """
}

process pathway_abund {
    echo true 
    publishDir 'picrust2_tmp', mode: 'copy'
    input:
    file pred_EC from output_EC 
    output:
    file "pathways_out/path_abun_strat.tsv.gz" into path_inferred
    script:
        """
        pathway_pipeline.py -i ${pred_EC}\
        -o pathways_out \
        -p ${params.threads}
        """
    stub:
        """
        mkdir -p pathways_out
        touch pathways_out/path_abun_strat.tsv.gz
        """
}

process add_description {
    echo true 
    publishDir 'picrust2_tmp', mode: 'copy'
    input: 
    file path_ab from path_inferred 
    file pred_EC from output_EC_desc 
    output: 
    file "pathways_out/path_abun_strat_desc.tsv.gz" 
    file "EC_metagenome_out/pred_metagenome_strat_desc.tsv.gz"
    script:
    """
    add_descriptions.py -i ${path_ab} \
    -m EC \
    -o pathways_out/path_abun_strat_desc.tsv.gz
    
    add_descriptions.py -i ${pred_EC} \
    -m METACYC \
    -o EC_metagenome_out/pred_metagenome_strat_desc.tsv.gz
    """
    stub:
    """
    mkdir -p pathways_out/
    mkdir -p EC_metagenome_out/
    touch pathways_out/path_abun_strat_desc.tsv.gz
    touch EC_metagenome_out/pred_metagenome_strat_desc.tsv.gz
    """
}

