rule compile_counts:
    input:
        expand("samples/htseq_count/{sample}_htseq_gene_count.txt", sample=SAMPLES),
    output:
        "results/{project_name}_mrna_counts.txt".format(project_name=config["project_name"]),
    run:
        import pandas as pd
        data = pd.DataFrame() 
        rownames =[]

        for file in input:
            temp_df=pd.read_csv(file, header=None, sep='\t')
            #append the rows of data
            data[file] = temp_df.iloc[:,1].values
            #append the rownames to double check
            rownames.append(temp_df.iloc[:,0].values)

        if all([(x==rownames[0]).all() for x in rownames]):
            data["gene_names"]= rownames[0]
            data.to_csv(output[0])
        else:
            raise Exception("rows are either not sorted or not identical please check") 

rule genecount:
    input:
        config["gscratch_path"] + "samples/genecounts_rmdp/{sample}_bam/{sample}_sort.rmd.bam"
    output:
        "samples/htseq_count/{sample}_htseq_gene_count.txt",
    params:
        gtf = config["gtf_file"]
    conda:
        "../envs/omic_qc_wf.yaml"
    #group: "short_post_align"
    shell:
        """htseq-count \
                -f bam \
                -r name \
                -s reverse \
                -m intersection-strict \
                {input} \
                {params.gtf} > {output}"""

rule sort:
    input:
        config["gscratch_path"] + "samples/genecounts_rmdp/{sample}_bam/{sample}.rmd.bam"
    output:
        temp(config["gscratch_path"] + "samples/genecounts_rmdp/{sample}_bam/{sample}_sort.rmd.bam")
    #group: "short_post_align"
    params:
        samtools="/opt/installed/samtools-1.6/bin/samtools"
    shell:
      """{params.samtools} sort -O bam -n {input} -o {output}"""

#todo switch to another duplicate caller?
rule picard:
    input:
        config["gscratch_path"] + "samples/star/{sample}_bam/Aligned.sortedByCoord.out.bam"
    output:
        temp(config["gscratch_path"] + "samples/genecounts_rmdp/{sample}_bam/{sample}.rmd.bam")
    #group: "short_post_align"
    params:
        picard=config["picard_tool"],
        metrics=config["gscratch_path"] + "samples/genecounts_rmdp/{sample}_bam/{sample}.rmd.metrics.text"
    resources:
        mem_mb=3500         
    #defined manually as 3 gigabytes in the rule below
    shell:
        """
        java -Xmx3g -jar {params.picard} \
        INPUT={input} \
        OUTPUT={output} \
        METRICS_FILE={params.metrics} \
        REMOVE_DUPLICATES=true
        """

rule index:
    input:
        config["gscratch_path"] + "samples/star/{sample}_bam/Aligned.sortedByCoord.out.bam"
    output:
        temp(config["gscratch_path"] + "samples/star/{sample}_bam/Aligned.sortedByCoord.out.bam.bai")
#    group: "short_post_align"
    params:
        samtools="/opt/installed/samtools-1.6/bin/samtools"
    shell:
        """{params.samtools} index {input} {output}"""


rule star:
    input:
        fwd = config["gscratch_path"] + "samples/trimmed/{sample}_R1_t.fq",
        rev = config["gscratch_path"] + "samples/trimmed/{sample}_R2_t.fq"
    output:
        temp(config["gscratch_path"] + "samples/star/{sample}_bam/Aligned.sortedByCoord.out.bam"),
        temp(config["gscratch_path"] + "samples/star/{sample}_bam/ReadsPerGene.out.tab"),
    threads: 5
    resources:
        mem_mb=55000
    params:
        pathToGenomeIndex = config["star_index"],
        gtf = config["gtf_file"],
        star = config["star_tool"],
        output_name_prefix = config["gscratch_path"] + "samples/star/{sample}_bam/",
    shell:
        """{params.star} --runThreadN {threads} --runMode alignReads --genomeDir {params.pathToGenomeIndex} \
        --readFilesIn {input.fwd} {input.rev} \
        --outFileNamePrefix {params.output_name_prefix} \
        --sjdbGTFfile {params.gtf} --quantMode GeneCounts \
        --sjdbGTFtagExonParentGene gene_name \
        --outSAMtype BAM SortedByCoordinate \
        """

rule fastqc:
    input:
        fwd = "samples/trimmed/{sample}_R1_t.fq",
        rev = "samples/trimmed/{sample}_R2_t.fq"
    output:
        fwd = "samples/fastqc/{sample}/{sample}_R1_t_fastqc.zip",
        rev = "samples/fastqc/{sample}/{sample}_R2_t_fastqc.zip"

    conda:
        "../envs/fastqc.yaml"
    shell:
        """fastqc --outdir samples/fastqc/{wildcards.sample} --extract  -f fastq {input.fwd} {input.rev}"""


rule trimming:
    input:
        fwd = "samples/raw/{sample}_R1.fq",
        rev = "samples/raw/{sample}_R2.fq"
    output:
        fwd = pipe(config["gscratch_path"] + "samples/trimmed/{sample}_R1_t.fq"),
        rev = pipe(config["gscratch_path"] + "samples/trimmed/{sample}_R2_t.fq"),
        single = temp("samples/trimmed/{sample}_R1_singletons.fq")
    params:
        sickle = config["sickle_tool"]
    shell:
        """{params.sickle} pe -f {input.fwd} -r {input.rev}  -l 40 -q 20 -t sanger  -o {output.fwd} -p {output.rev} -s {output.single} &> {input.fwd}.log
        """

