rule Bowtie2_Index:
    output:
        expand(
            f"{config['bowtie2_patho_db']}{{ext}}",
            ext=[
                ".1.bt2l",
                ".2.bt2l",
                ".3.bt2l",
                ".4.bt2l",
                ".rev.1.bt2l",
                ".rev.2.bt2l",
            ],
        ),
    input:
        ref=config["bowtie2_patho_db"],
    conda:
        "../envs/bowtie2.yaml"
    envmodules:
        *config["envmodules"]["bowtie2"],
    threads: 1
    log:
        f"logs/BOWTIE2_BUILD/{config['bowtie2_patho_db']}.log",
    shell:
        "bowtie2-build-l {input.ref} {input.ref} > {log} 2>&1"


rule Bowtie2_Pathogenome_Alignment:
    output:
        bam="results/BOWTIE2/{sample}/AlignedToPathogenome.bam",
        bai="results/BOWTIE2/{sample}/AlignedToPathogenome.bam.bai",
    input:
        fastq="results/CUTADAPT_ADAPTER_TRIMMING/{sample}.trimmed.fastq.gz",
        db=rules.Bowtie2_Index.output,
    params:
        PATHO_DB=lambda wildcards, input: config["bowtie2_patho_db"],
    threads: 10
    log:
        "logs/BOWTIE2/{sample}.log",
    conda:
        "../envs/bowtie2.yaml"
    envmodules:
        *config["envmodules"]["bowtie2"],
    benchmark:
        "benchmarks/BOWTIE2/{sample}.benchmark.txt"
    message:
        "Bowtie2_Pathogenome_Alignment: ALIGNING SAMPLE {input.fastq} TO PATHOGENOME WITH BOWTIE2"
    shell:
        """bowtie2 --large-index -x {params.PATHO_DB} --end-to-end --threads {threads} --very-sensitive -U {input.fastq} 2> {log} | samtools view -bS -q 1 -h -@ {threads} - | samtools sort - -@ {threads} -o {output.bam} >> {log};"""
        """samtools index {output.bam}"""
