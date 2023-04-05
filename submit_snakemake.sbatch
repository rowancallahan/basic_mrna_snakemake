#!/bin/bash
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=3gb
#SBATCH --time=17:00:00
#SBATCH --partition exacloud  
#SBATCH --job-name workflow_submission

snakemake --configfile config.yaml --profile profile --keep-going --rerun-incomplete

