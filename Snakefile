from snakemake.utils import min_version
min_version("6.0")

configfile: "config/config.yaml"

module dada2:
    snakefile: "workflow/rules/dada2.smk"
    config: config['dada2']
module picrust:
    snakefile: "workflow/rules/picrust.smk"
    

use rule * from dada2 as dada2_*
use rule * from picrust as picrust_*