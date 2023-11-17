#!/bin/bash
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=3gb
#SBATCH --time=36:00:00
#SBATCH --partition exacloud  
#SBATCH --job-name workflow_submission
#SBATCH --output=logs/workflow_submission_%j.log
#SBATCH --mail-user=spilioto@ohsu.edu
#SBATCH --mail-type=FAIL

snakemake --configfile config.yaml --profile profile --keep-going --rerun-incomplete

