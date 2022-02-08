rule FastQC_BeforeTrimming:
    """Run fastq before trimming"""
    output:
        html="results/FASTQC_BEFORE_TRIMMING/{sample}_fastqc.html",
        zip="results/FASTQC_BEFORE_TRIMMING/{sample}_fastqc.zip",
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
        zip="results/FASTQC_AFTER_TRIMMING/{sample}.trimmed_fastqc.zip",
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


rule MultiQC:
    """Run MultiQC"""
    output:
        html="results/MULTIQC/multiqc_report.html",
    input:
        unpack(multiqc_input),
    log:
        "logs/MULTIQC/MULTIQC.log",
    conda:
        "../envs/multiqc.yaml"
    params:
        config=os.path.join(WORKFLOW_DIR, "envs", "multiqc_config.yaml"),
    envmodules:
        *config["envmodules"]["MultiQC"],
    benchmark:
        "benchmarks/MULTIQC/MULTIQC.benchmark.txt"
    message:
        "COMBINING QUALITY CONTROL METRICS WITH MULTIQC"
    shell:
        'echo {input} | tr " " "\n" > {output.html}.fof;'
        "multiqc -c {params.config} -l {output.html}.fof --verbose --force --outdir results/MULTIQC &> {log}"
