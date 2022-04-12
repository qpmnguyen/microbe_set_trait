from snakemake.utils import min_version
min_version("6.0")

module picrust:
    snakefile: "workflow/rules/picrust.smk"

module db_prep:
    snakefile: "workflow/rules/db_prep.smk"
    config: config

use rule * from picrust as picrust_*
use rule * from db_prep as db_prep_