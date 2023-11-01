rule spliceq:
    input:
        config["gscratch_path"] + "samples/genecounts_rmdp/{sample}_bam/{sample}.sortposition.rmd.bam"
    output:
        "samples/spliceq/{sample}_counts.tsv"
    params:
        gtf=config["gtf_file"],
    resources:
        mem_mb=32000
    threads: 4 
    shell:
        """
        SPLICE-q.py --NrProcesses 4 --gtffile {params.gtf} --bamfile {input} --outfile {output}
        """

rule sort_position:
    input:
        config["gscratch_path"] + "samples/genecounts_rmdp/{sample}_bam/{sample}_sort.rmd.bam"
    output:
        temp(config["gscratch_path"] + "samples/genecounts_rmdp/{sample}_bam/{sample}.sortposition.rmd.bam")
    #group: "short_post_align"
    params:
        samtools="/opt/installed/samtools-1.6/bin/samtools"
    shell:
      """{params.samtools} sort -O bam {input} -o {output}"""



