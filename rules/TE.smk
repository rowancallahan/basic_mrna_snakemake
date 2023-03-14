rule combine_TE:
    input:
        expand("samples/te_count/{sample}_path_to_out_file", sample=SAMPLES)
    output:
        expand("results/{project_name}_TE_counts.txt", project_name=project_title),
    run:
        """
        python code to combine all of the different files together
        """

#todo figoure out what the actual otputs are for outfile
rule TEtranscripts:
    input:
        "samples/star_TE/{sample}/Aligned.out.bam",
    output:
        "samples/te_count/{sample}_path_to_out_file",
    conda:
        "../envs/TE.yaml"
    params:
        gtf=config["gtf_file"],
        TE_gtf=config["TE_gtf"],
        outdir="samples/te_count",
    shell:
        """
        TEcount --format BAM --GTF {params.gtf} --TE_gtf {params.TE_gtf} \
        --mode multi -b {input} --project {wildcards.sample} \
        --outdir {params.outdir}
        """

rule map_TE:
    input:
        fwd = "samples/trimmed/{sample}_R1_t.fq",
        rev = "samples/trimmed/{sample}_R2_t.fq"
    output:
        "samples/star_TE/{sample}/Aligned.out.bam",
        "samples/star_TE/{sample}/Log.final.out"
    params:
        index=config["star_index"],
        gtf=config["gtf_file"]
    run:
        STAR=config["star_tool"]

        shell("""
                {STAR} --runThreadN 12 --genomeDir {params.index} --sjdbGTFfile {params.gtf} \
                       --sjdbOverhang 100 --readFilesIn {input.fwd} {input.rev} \
                       --outFileNamePrefix samples/star_TE/{wildcards.sample}/ \
                       --outSAMtype BAM Unsorted --winAnchorMultimapNmax 200 --outFilterMultimapNmax 100""")


