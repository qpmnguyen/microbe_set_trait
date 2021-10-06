#!/usr/bin/env nextflow

biom_file = Channel.fromPath("data/gevers_biom.biom")
seqs_file = Channel.fromPath("gevers_seq.fa")

process place_seqs {
    echo true
    publishDir 'picrust2_tmp/'
    input:
    file seqs from seqs_file
    output: 
    file "picrust2_tmp/placed_seqs.tre" into placed_tree
    """
    place_seqs.py -s $seqs -o placed_seqs.tre --threads 2 -a hmmalign --print_cmds
    """
}

process hidden_state {
    echo true
    publishDir 'picrust2_tmp/'
    input: 
    file tree from placed_tree
    output:
    file "picrust2_tmp/pred_traits_marker.tsv.gz" into pred_marker 
    file "picrust2_tmp/pred_traits_EC.tsv.gz" into pred_EC 
    """
    hsp.py -t $tree -i 16S -o picrust2_tmp/pred_traits_marker -m mp --check -p 2 -n
    hsp.py -t $tree -i EC -o picrust2_tmp/pred_traits_EC -m mp --check -p 2
    """
}

process metagene_pred {
    echo true 
    publishDir 'picrust2_tmp'
    input: 
    file biom from biom_file 
    file marker from pred_marker
    file enzyme from pred_EC 
    output: 
    file 'picrust2_tmp/EC_metagenome_out/*' into output_EC
    """
    metagenome_pipeline.py -i \
    -m ${pred_marker} \
    -f ${pred_EC} \
    -o picrust2_tmp/EC_metagenome_out \
    --strat_out \
    """
}

process pathway_abund {
    """
    pathway_pipeline.py -i \
    -o \
    """
}

process add_description {

}

