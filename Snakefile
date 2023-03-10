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
configfile:"config.yaml"
project_title= config["project_name"]


##############################
# FILEPATH LINKING HERE \/ \/#
##############################
SAMPLES, = glob_wildcards("samples/raw/{sample}_R1.fq")
print(SAMPLES)



rule all:
    input:
        #expand("./results/{project_name}_mrna_counts.txt", project_name=project_title),
        expand("results/{project_name}_TE_counts.txt", project_name=project_title),
        #expand("./results/{project_name}_quality_control_metrics.txt", project_name=project_title),
        #expand("./results/{project_name}_circ_counts.txt", project_name=project_title),
        #expand("./results/{project_name}_microbe_counts.txt", project_name=project_title),


        

include: "rules/mrna.smk"
include: "rules/qc.smk"
include: "rules/TE.smk"
include: "rules/circ.smk"
include: "rules/microbe.smk"
