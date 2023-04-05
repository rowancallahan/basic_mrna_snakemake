rule combine_TE:
    input:
        unstranded = expand("samples/te_count/{sample}.cntTable", sample=SAMPLES),
        rev = expand("samples/te_count_rev/{sample}.cntTable", sample=SAMPLES),
    output:
        unstranded = "results/{project_name}_unstranded_TE_counts.txt".format(project_name=config["project_name"]),
        rev =  "results/{project_name}_reverse_TE_counts.txt".format(project_name=config["project_name"]),

    run:
        import pandas as pd

        unstranded= pd.DataFrame() 
        rownames =[]
        for file in input.unstranded:
            temp_df=pd.read_csv(file, header=0, sep='\t')
            unstranded[file] = temp_df.iloc[:,1].values
            if not rownames:
                unstranded["gene_name"] =  temp_df.iloc[:,0].values
        unstranded.to_csv(output.unstranded)
        
        rownames= []
        rev= pd.DataFrame() 
        for file in input.rev:
            temp_df=pd.read_csv(file, header=0, sep='\t')
            rev[file] = temp_df.iloc[:,1].values
            if not rownames:
                rev["gene_name"] =  temp_df.iloc[:,0].values
        rev.to_csv(output.rev)


rule TEtranscripts_rev:
    input:
        config["gscratch_path"] + "samples/star_TE/{sample}/Aligned.sortedByCoord.out.bam",
    output:
        "samples/te_count_rev/{sample}.cntTable"
    params:
        gtf=config["gtf_file"],
        TE_gtf=config["TE_gtf"],
        outdir="samples/te_count_rev",
    resources:
        mem_mb=32000
    shell:
        """
        TEcount --format BAM --GTF {params.gtf} --TE {params.TE_gtf} \
        --mode multi -b {input} --project {wildcards.sample} \
        --outdir {params.outdir} --sortByPos --stranded reverse
        """

rule TEtranscripts:
    input:
        config["gscratch_path"] + "samples/star_TE/{sample}/Aligned.sortedByCoord.out.bam",
    output:
        "samples/te_count/{sample}.cntTable"
    params:
        gtf=config["gtf_file"],
        TE_gtf=config["TE_gtf"],
        outdir="samples/te_count",
    resources:
        mem_mb=32000
    shell:
        """
        TEcount --format BAM --GTF {params.gtf} --TE {params.TE_gtf} \
        --mode multi -b {input} --project {wildcards.sample} \
        --outdir {params.outdir} --sortByPos
        """

rule map_TE:
    input:
        fwd = config["gscratch_path"] + "samples/clumped/{sample}_R1_clumped.fq",
        rev = config["gscratch_path"] + "samples/clumped/{sample}_R2_clumped.fq"
    output:
        config["gscratch_path"] + "samples/star_TE/{sample}/Aligned.sortedByCoord.out.bam",
        config["gscratch_path"] + "samples/star_TE/{sample}/Log.final.out",
    params:
        index=config["star_index"],
        gtf=config["gtf_file"],
        outfile_prefix = config["gscratch_path"] + "samples/star_TE/{sample}/",
        STAR=config["star_tool"],

    threads: 8 
    resources:
        mem_mb=55000
    shell:
        """
        {params.STAR} --runThreadN {threads} \
        --genomeDir {params.index} \
        --sjdbGTFfile {params.gtf} \
        --readFilesIn {input.fwd} {input.rev} \
        --outFileNamePrefix {params.outfile_prefix} \
        --outSAMtype BAM SortedByCoordinate \
        --winAnchorMultimapNmax 100 \
        --outFilterMultimapNmax 100
        """


rule clumpify:
    input:
        fwd = config["gscratch_path"] + "samples/trimmed/{sample}_R1_t.fq",
        rev = config["gscratch_path"] + "samples/trimmed/{sample}_R2_t.fq",
    output:
        fwd = temp(config["gscratch_path"] + "samples/clumped/{sample}_R1_clumped.fq"),
        rev = temp(config["gscratch_path"] + "samples/clumped/{sample}_R2_clumped.fq"),
    resources:
        mem_mb=66000
    shell:
        """
        export PATH=$PATH:/home/groups/CEDAR/callahro/software/bbmap/bbmap/
        clumpify.sh in={input.fwd} in2={input.rev} dedupe=t out={output.fwd} out2={output.rev}
        """


