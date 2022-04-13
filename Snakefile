from snakemake.utils import min_version
min_version("6.0")

module dada2:
    snakefile: "workflow/rules/dada2.smk"
module picrust:
    snakefile: "workflow/rules/picrust.smk"

use rule * from dada2 as dada2_*
use rule * from picrust as picrust_*