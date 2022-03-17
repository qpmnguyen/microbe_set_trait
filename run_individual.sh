#!/usr/bin/env bash

# Quang Nguyen - 03/14/2022
# This script is used interactively to run individual scripts via nohup
source /optnfs/common/miniconda3/etc/profile.d/conda.sh
conda activate microbe_trait 

Rscript -e "workflowr::wflow_publish('analysis/agp_download.Rmd', verbose=TRUE)" 
