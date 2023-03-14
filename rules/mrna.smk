rule compile_counts:
    input:
        expand("samples/htseq_count/{sample}_htseq_gene_count.txt",sample=SAMPLES)
    output:
        expand("results/{project_name}_mrna_counts.txt", project_name=project_title),
    script:
        "../scripts/compile_counts_table.py"


rule genecount:
    input:
        "samples/genecounts_rmdp/{sample}_bam/{sample}_sort.rmd.bam"
    output:
        "samples/htseq_count/{sample}_htseq_gene_count.txt",

    log:
        "logs/genecount/{sample}_genecount.log"
    params:
        name = "genecount_{sample}",
        gtf = config["gtf_file"]
    conda:
        "../envs/omic_qc_wf.yaml"
    threads: 1
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
      "samples/genecounts_rmdp/{sample}_bam/{sample}.rmd.bam"
    output:
      "samples/genecounts_rmdp/{sample}_bam/{sample}_sort.rmd.bam"
    params:
      name = "sort_{sample}",
      mem = "6400"

    conda:
      "../envs/omic_qc_wf.yaml"
    shell:
      """samtools sort -O bam -n {input} -o {output}"""

rule picard:
    input:
        "samples/star/{sample}_bam/Aligned.sortedByCoord.out.bam"
    output:
        temp("samples/genecounts_rmdp/{sample}_bam/{sample}.rmd.bam")
    params:
        name="rmd_{sample}",
        mem="5300"
    run:
      picard=config["picard_tool"]

      shell("java -Xmx3g -jar {picard} \
      INPUT={input} \
      OUTPUT={output} \
      METRICS_FILE=samples/genecounts_rmdp/{wildcards.sample}_bam/{wildcards.sample}.rmd.metrics.text \
      REMOVE_DUPLICATES=true")

rule index:
    input:
        "samples/star/{sample}_bam/Aligned.sortedByCoord.out.bam"
    output:
        "samples/star/{sample}_bam/Aligned.sortedByCoord.out.bam.bai"

    conda:
        "../envs/omic_qc_wf.yaml"
    shell:
        """samtools index {input} {output}"""


rule star:
    input:
        fwd = "samples/trimmed/{sample}_R1_t.fq",
        rev = "samples/trimmed/{sample}_R2_t.fq"
    output:
        "samples/star/{sample}_bam/Aligned.sortedByCoord.out.bam",
        "samples/star/{sample}_bam/ReadsPerGene.out.tab",
        "samples/star/{sample}_bam/Log.final.out"
    params:
        gtf=config["gtf_file"]

    run:
         STAR=config["star_tool"],
         pathToGenomeIndex = config["star_index"]

         shell("""
                {STAR} --runThreadN {threads} --runMode alignReads --genomeDir {pathToGenomeIndex} \
                --readFilesIn {input.fwd} {input.rev} \
                --outFileNamePrefix samples/star/{wildcards.sample}_bam/ \
                --sjdbGTFfile {params.gtf} --quantMode GeneCounts \
                --sjdbGTFtagExonParentGene gene_name \
                --outSAMtype BAM SortedByCoordinate \
                """)

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
        fwd = "samples/trimmed/{sample}_R1_t.fq",
        rev = "samples/trimmed/{sample}_R2_t.fq",
        single = "samples/trimmed/{sample}_R1_singletons.fq"
    run:
        sickle = config["sickle_tool"]
        shell("{sickle} pe -f {input.fwd} -r {input.rev}  -l 40 -q 20 -t sanger  -o {output.fwd} -p {output.rev} -s {output.single} &> {input.fwd}.log")

