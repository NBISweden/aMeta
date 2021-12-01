rule FastQC_BeforeTrimming:
    """Run fastq before trimming"""
    output:
        html="results/FASTQC_BEFORE_TRIMMING/{sample}_fastqc.html",
    input:
        fastq=lambda wildcards: samples.loc[wildcards.sample].fastq,
    conda:
        "../envs/fastqc.yaml"
    envmodules:
        *config["envmodules"]["FastQC_BeforeTrimming"],
    benchmark:
        "benchmarks/FASTQC_BEFORE_TRIMMING/{sample}.benchmark.txt"
    message:
        "RUNNING QUALITY CONTROL WITH FASTQC FOR SAMPLE {input.fastq} BEFORE TRIMMING ADAPTERS"
    log:
        "logs/FASTQC_BEFORE_TRIMMING/{sample}.log",
    threads: 1
    shell:
        "fastqc {input.fastq} --threads {threads} --nogroup --outdir results/FASTQC_BEFORE_TRIMMING &> {log}"


rule MultiQC_BeforeTrimming:
    output:
        html="results/MULTIQC_BEFORE_TRIMMING/multiqc_report.html",
    input:
        expand("results/FASTQC_BEFORE_TRIMMING/{sample}_fastqc.html", sample=SAMPLES),
    log:
        "logs/MULTIQC_BEFORE_TRIMMING/MULTIQC_BEFORE_TRIMMING.log",
    conda:
        "../envs/multiqc.yaml"
    envmodules:
        *config["envmodules"]["MultiQC_BeforeTrimming"],
    benchmark:
        "benchmarks/MULTIQC_BEFORE_TRIMMING/MULTIQC_BEFORE_TRIMMING.benchmark.txt"
    message:
        "COMBINING QUALITY CONTROL METRICS WITH MULTIQC BEFORE TRIMMING ADAPTERS"
    shell:
        "multiqc results/FASTQC_BEFORE_TRIMMING --verbose --force --outdir results/MULTIQC_BEFORE_TRIMMING &> {log}"


rule Cutadapt_Adapter_Trimming:
    output:
        fastq="results/CUTADAPT_ADAPTER_TRIMMING/{sample}.trimmed.fastq.gz",
    input:
        fastq=lambda wildcards: samples.loc[wildcards.sample].fastq,
    params:
        illumina_universal_adapter=config["illumina_universal_adapter"],
    log:
        "logs/CUTADAPT_ADAPTER_TRIMMING/{sample}.log",
    conda:
        "../envs/cutadapt.yaml"
    envmodules:
        *config["envmodules"]["Cutadapt_Adapter_Trimming"],
    benchmark:
        "benchmarks/CUTADAPT_ADAPTER_TRIMMING/{sample}.benchmark.txt"
    message:
        "TRIMMING ADAPTERS FOR SAMPLE {input.fastq} WITH CUTADAPT"
    threads: 1
    shell:
        "cutadapt -a {params.illumina_universal_adapter} --minimum-length 30 -o {output.fastq} {input.fastq} &> {log}"


rule FastQC_AfterTrimming:
    output:
        html="results/FASTQC_AFTER_TRIMMING/{sample}.trimmed_fastqc.html",
    input:
        fastq="results/CUTADAPT_ADAPTER_TRIMMING/{sample}.trimmed.fastq.gz",
    threads: 1
    log:
        "logs/FASTQC_AFTER_TRIMMING/{sample}.log",
    conda:
        "../envs/fastqc.yaml"
    envmodules:
        *config["envmodules"]["FastQC_AfterTrimming"],
    benchmark:
        "benchmarks/FASTQC_AFTER_TRIMMING/{sample}.benchmark.txt"
    message:
        "RUNNING QUALITY CONTROL WITH FASTQC FOR SAMPLE {input.fastq} AFTER TRIMMING ADAPTERS"
    shell:
        "fastqc {input.fastq} --threads {threads} --nogroup --outdir results/FASTQC_AFTER_TRIMMING &> {log}"


# FIXME? Also add cutadapt output
rule MultiQC_AfterTrimming:
    """Run MultiQC on trimmed data"""
    output:
        html="results/MULTIQC_AFTER_TRIMMING/multiqc_report.html",
    input:
        expand(
            "results/FASTQC_AFTER_TRIMMING/{sample}.trimmed_fastqc.html",
            sample=SAMPLES,
        ),
    log:
        "logs/MULTIQC_AFTER_TRIMMING/MULTIQC_AFTER_TRIMMING.log",
    conda:
        "../envs/multiqc.yaml"
    envmodules:
        *config["envmodules"]["MultiQC_AfterTrimming"],
    benchmark:
        "benchmarks/MULTIQC_AFTER_TRIMMING/MULTIQC_AFTER_TRIMMING.benchmark.txt"
    message:
        "COMBINING QUALITY CONTROL METRICS WITH MULTIQC AFTER TRIMMING ADAPTERS"
    shell:
        "multiqc results/FASTQC_AFTER_TRIMMING --verbose --force --outdir results/MULTIQC_AFTER_TRIMMING &> {log}"
