__author__ = "Rowan Callahan"
__email__ = "callahro@ohsu.edu"
__license__ = "MIT"

#data should have raw data and intermediate tabulations  
#results should have tables and main outputs
#all graphing will mainly be done using external scripts to reduce dependencies
#everything should exist inside the main rules specified below

##############################
# CONFIG FILE LOADING HERE \/#
##############################
configfile:
  "config.yaml"

##############################
# FILEPATH LINKING HERE \/ \/#
##############################
SAMPLES, = glob_wildcards("samples/raw/{sample}_R1.fq")
print(SAMPLES)

rule all:
    input:
        "results/{project_name}_mrna_counts.txt".format(project_name=config["project_name"]),

include: "rules/mrna.smk"


