rule telescope:
    input:
        config["gscratch_path"] + "samples/star_TE/{sample}/Aligned.sortedByCoord.out.bam",
    output:
        "samples/telescope_count/{sample}/telescope-TE_counts.tsv""
    params:
        TE_gtf=config["TE_gtf"],
        outdir="samples/telescope_count/{sample}",
        sorted_bam="samples/telescope_count/{sample}.bam",
    resources:
        mem_mb=32000
    log:
        "logs/telescope/{sample}_telescope.log"
    threads: 4 
    shell:
        """
        samtools sort -n {input} -o {params.sorted_bam}
        telescope assign --ncpu {threads} --outdir {params.outdir} {params.sorted_bam} {params.TE_gtf}
        rm {params.sorted_bam}
        """


