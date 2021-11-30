rule Bowtie2_Pathogenome_Alignment:
    input:
        fastq = "results/CUTADAPT_ADAPTER_TRIMMING/{sample}.trimmed.fastq.gz"
    params:
        # PATHO_DB = "/proj/snic2018-8-150/uppstore2018095/private/NBIS_Demo/PathoGenome/library.pathogen.fna"
        PATHO_DB=config["bowtie2_patho_db"]
    threads: 10
    log:
        "logs/BOWTIE2/{sample}.log"
    conda:
        "../envs/bowtie2.yaml"
    envmodules:
        *config["envmodules"]["Bowtie2_Pathogenome_Alignment"],
    benchmark:
        "benchmarks/BOWTIE2/{sample}.benchmark.txt"
    message:
        "ALIGNING SAMPLE {input.fastq} TO PATHOGENOME WITH BOWTIE2"
    output:
        bam = "results/BOWTIE2/{sample}/AlignedToPathogenome.bam"
    run:
        "bowtie2 --large-index -x {params.PATHO_DB} --end-to-end --threads {threads} --very-sensitive -U {input.fastq} | samtools view -bS -q 1 -h -@ {threads} - | samtools sort - -@ {threads} > {output.bam} &> {log}; "
        "samtools index {output.bam}"
