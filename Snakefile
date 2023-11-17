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

##############################
# FILEPATH LINKING HERE \/ \/#
##############################
SAMPLES, = glob_wildcards("samples/raw/{sample}_R1.fq")
print(SAMPLES)
gscratch_location = "/home/exacloud/gscratch/CEDAR/cfRNA-seq-pipeline/multi_omic_pipeline/"

rule all:
    input:
        #"results/{project_name}_mrna_counts.txt".format(project_name=config["project_name"]),
        #"results/{project_name}_unstranded_TE_counts.txt".format(project_name=config["project_name"]),
        #expand("./results/{project_name}_quality_control_metrics.txt", project_name=config["project_name"]),
        expand("./results/{project_name}_circ_CE.txt", project_name=config["project_name"]),
        expand("./results/{project_name}_circ_CQ.txt", project_name=config["project_name"]),
        #expand("./results/{project_name}_microbe_counts.txt", project_name=config["project_name"]),
        #expand("samples/spliceq/{sample}_counts.tsv", sample=SAMPLES),


include: "rules/mrna.smk"
##include: "rules/TE.smk"
##include: "rules/spliceQ.smk"
include: "rules/circ.smk"
##include: "rules/microbe.smk"
##include: "rules/qc.smk"
