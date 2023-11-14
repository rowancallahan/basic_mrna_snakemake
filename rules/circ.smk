rule ciriquant:
        input:
            R1 = config["gscratch_path"] + "samples/trimmed/{sample}_R1_t.fq",
            R2 = config["gscratch_path"] + "samples/trimmed/{sample}_R2_t.fq"
        output:
            config["gscratch_path"] + "samples/ciriquant/{sample}/{sample}.gtf",
            temp(config["gscratch_path"] + "samples/ciriquant/{sample}/align/{sample}.sorted.bam"),
            temp(config["gscratch_path"] + "samples/ciriquant/{sample}/circ/{sample}_denovo.sorted.bam")
        params:
            ciri_conf = config["ucsc_ciri_conf"]
            output_name_prefix = config["gscratch_path"] + "samples/ciriquant/{sample}/"
        conda:
            "../envs/ciriquant.yaml"
        threads: 8
        resources:
            mem_mb=64000
        shell:
            """CIRIquant -t {threads} -1 {input.R1} -2 {input.R2} --config {params.ciri_conf} -o {params.output_name_prefix} -p {wildcards.sample} """

rule star_CE:
    input:
        fwd = config["gscratch_path"] + "samples/trimmed/{sample}_R1_t.fq",
        rev = config["gscratch_path"] + "samples/trimmed/{sample}_R2_t.fq"
    output:
        temp(config["gscratch_path"] + "samples/star_exp/{sample}_bam/Aligned.sortedByCoord.out.bam"),
        config["gscratch_path"] + "samples/star_exp/{sample}_bam/ReadsPerGene.out.tab",
        config["gscratch_path"] + "samples/star_exp/{sample}_bam/Chimeric.out.junction"
    params:
        ref_gtf_exp=config["ref_gtf_exp"],
        STAR=config["star_tool"],
        pathToGenomeIndex = config["star_index_exp"],
        output_name_prefix = config["gscratch_path"] + "samples/star_exp/{sample}_bam/"
    threads: 8
    resources:
        mem_mb=64000
    shell:
        """{params.STAR} --runThreadN {threads} --runMode alignReads --genomeDir {pathToGenomeIndex} \
                --readFilesIn {input.fwd} {input.rev} \
                --outFileNamePrefix {params.output_name_prefix} \
                --sjdbGTFfile {params.ref_gtf_exp} --quantMode GeneCounts \
                --sjdbGTFtagExonParentGene gene_name \
                --outSAMtype BAM SortedByCoordinate \
                --chimSegmentMin 10 \
                #--readFilesCommand zcat \
                --twopassMode Basic
                """)


rule circ_explorer_parse_star:
        input: config["gscratch_path"] + "samples/star_exp/{sample}_bam/Chimeric.out.junction"
        output: config["gscratch_path"] + "samples/circ_explorer_star/{sample}/{sample}_backspliced.junction.bed"
        threads: 4
        resources:
            mem_mb=10000
        conda:
                "../envs/circexplorer.yaml"
        shell: """CIRCexplorer2 parse -t STAR {input} -b {output}"""

rule circ_explorer_annotate_star:
        input: config["gscratch_path"] + "samples/circ_explorer_star/{sample}/{sample}_backspliced.junction.bed"
        output: config["gscratch_path"] + "samples/circ_explorer_star/{sample}_known_circ.txt"
        params:
                ref_fa_exp = config["ref_fa_exp"],
                ref_gtf_exp = config["ref_txt_exp"]
        conda:
                "../envs/circexplorer.yaml"
        threads: 2
        resources:
            mem_mb=10000
        shell:
            """ CIRCexplorer2 annotate -r {params.ref_gtf_exp} -g {params.ref_fa_exp} -b {input} -o {output} """

rule circ_compiler_CE:
        input:
            expand(config["gscratch_path"] + "samples/circ_explorer_star/{sample}_known_circ.txt", sample = SAMPLES)
        output:
            expand("./results/{project_name}_circ_CE.txt", project_name=project_title)
        script:
            "./scripts/compile_CE.R"

rule circ_compiler_CQ:
        input:
            expand(config["gscratch_path"] + "samples/ciriquant/{sample}/{sample}.gtf", sample = SAMPLES)
        output:
            expand("./results/{project_name}_circ_CQ.txt", project_name=project_title)
        script:
            "./scripts/compile_CQ.R"
