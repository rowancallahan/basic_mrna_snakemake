## MultiOMIC analysis of cfRNA sequencing pipeline

Rowan Callahan and Elias Spiliotopolous

THe following pipeline has around four different parts:

# 1: Snakemake 
 - This will have our main snakefile and our five different sub rules that we use    
# 2: Environments
 - all of the config yaml files that we use w conda to setup environments
# 3: Output Folders
 - data for all of the intermediary results
 - results for all of our tabulated results 
# 4: Config files
 - config.yaml here
 - config.yaml in the profile
 - metadata.txt metadata for all of the samples
 - 


TODO:

1: get samples names from the metadata not from the file locations need to adapt the file linking script that we are currently using
2: Start minimum viable submittable job for pipeline overall
3: finish readme and config setup for the pipeline
4: setup for new pipeline usage script
