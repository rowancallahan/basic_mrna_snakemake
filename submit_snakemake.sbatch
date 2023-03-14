#!/bin/bash
#SBATCH --cpus-per-task=8 
#SBATCH --mem-per-cpu=3gb
#SBATCH --time=17:00:00
#SBATCH --partition exacloud  
#SBATCH --job-name workflow_submission

#snakemake -j 100 --use-conda --rerun-incomplete --restart-times 3 --latency-wait 20 --cluster-config cluster.json --cluster "sbatch -p {cluster.partition} -N {cluster.N}  -t {cluster.t} -o {cluster.o} -e {cluster.e} -J {cluster.J} -c {cluster.c} --mem {cluster.mem}" -s Snakefile

snakemake --configfile config.yaml --rerun-incomplete --use-conda --cluster "sbatch --cpus-per-task {threads} --mem {resources.mem_mb}m" --jobs 5 --latency-wait 20 

