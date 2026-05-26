rule Bowtie2_Index:
    output:
        expand(
            f"{config['bowtie2_db']}{{ext}}",
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
        ref=ancient(config["bowtie2_db"]),
    conda:
        "../envs/bowtie2.yaml"
    envmodules:
        *config["envmodules"]["bowtie2"],
    threads: 1
    log:
        f"{config['bowtie2_db']}_BOWTIE2_BUILD.log",
    shell:
        "bowtie2-build --large-index --threads {threads} {input.ref} {input.ref} > {log} 2>&1"


rule Bowtie2_Alignment:
    output:
        bam="results/BOWTIE2/{sample}/AlignedToBowtie2DB.bam",
        bai="results/BOWTIE2/{sample}/AlignedToBowtie2DB.bam.bai",
    input:
        fastq="results/CUTADAPT_ADAPTER_TRIMMING/{sample}.trimmed.fastq.gz",
        db=rules.Bowtie2_Index.output,
    params:
        BOWTIE2_DB=lambda wildcards, input: config["bowtie2_db"],
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
        "Bowtie2_Alignment: ALIGNING SAMPLE {input.fastq} WITH BOWTIE2"
    shell:
        """bowtie2 --large-index -x {params.BOWTIE2_DB} --end-to-end --threads {threads} --very-sensitive -U {input.fastq} > $(dirname {output.bam})/AlignedToBowtie2DB.sam 2> {log};"""
        """grep @ $(dirname {output.bam})/AlignedToBowtie2DB.sam | awk '!seen[$2]++' > $(dirname {output.bam})/header_nodups.txt;"""
        """grep -v '^@' $(dirname {output.bam})/AlignedToBowtie2DB.sam > $(dirname {output.bam})/AlignedToBowtie2DB.noheader.sam;"""
        """cat $(dirname {output.bam})/header_nodups.txt $(dirname {output.bam})/AlignedToBowtie2DB.noheader.sam > $(dirname {output.bam})/AlignedToBowtie2DB.nodups.sam;"""
        """samtools view -bS -q 1 -h -@ {threads} $(dirname {output.bam})/AlignedToBowtie2DB.nodups.sam | samtools sort - -@ {threads} -o {output.bam};"""
        """samtools index {output.bam};"""
        """rm $(dirname {output.bam})/header_nodups.txt;"""
        """rm $(dirname {output.bam})/AlignedToBowtie2DB.noheader.sam;"""
        """rm $(dirname {output.bam})/AlignedToBowtie2DB.nodups.sam;"""
        """rm $(dirname {output.bam})/AlignedToBowtie2DB.sam;"""
