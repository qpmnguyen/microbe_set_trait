from snakemake.utils import min_version
min_version("6.0")

configfile: "config/config.yaml"

module dada2:
    snakefile: "workflow/rules/dada2.smk"
    config: config['dada2']
module picrust:
    snakefile: "workflow/rules/picrust.smk"
module db2set:
    snakefile: "workflow/rules/db2set.smk"
    config: config 
module pred:
    snakefile: "workflow/rules/pred.smk"
    config: config['predict']
module trait_scores:
    snakefile: "workflow/rules/trait_scores.smk"
    config: config
module coverage:
    snakefile: "workflow/rules/coverage.smk"
    config: config['coverage']
module feat_importance:
    snakefile: "workflow/rules/feat_importance.smk"
    config: config['predict']

use rule * from dada2 as dada2_*
use rule * from picrust as picrust_*
use rule * from db2set as db2set_*
use rule * from pred as pred_*
use rule * from trait_scores as trait_scores_*
use rule * from coverage as coverage_*
use rule * from feat_importance as feat_importance_*