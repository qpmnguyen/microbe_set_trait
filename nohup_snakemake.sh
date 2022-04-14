#!/usr/bin/env bash

CONDA="conda activate notebooks"
eval $conda 

while getopts p:c: flag
do
    case "${flag}" in
        p) pipeline=${OPTARG};;
        c) cores=${OPTARG};;
    esac
done


if [ "$pipeline" = "dada2" ]; then
    target="dada2_all"
else
    echo "There is no $pipeline pipeline available"
    exit 1
fi

cmd="snakemake --use-conda --cores $cores $target"
echo "$cmd"
eval $cmd