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
        *config["envmodules"]["fastqc"],
    benchmark:
        "benchmarks/FASTQC_BEFORE_TRIMMING/{sample}.benchmark.txt"
    resources:
        mem_mb=10240,
    message:
        "FastQC_BeforeTrimming: RUNNING QUALITY CONTROL WITH FASTQC FOR SAMPLE {input.fastq} BEFORE TRIMMING ADAPTERS"
    log:
        "logs/FASTQC_BEFORE_TRIMMING/{sample}.log",
    threads: 2
    shell:
        "fastqc {input.fastq} --memory {resources.mem_mb} --threads {threads} --nogroup --outdir results/FASTQC_BEFORE_TRIMMING &> {log}"


rule Cutadapt_Adapter_Trimming:
    output:
        fastq="results/CUTADAPT_ADAPTER_TRIMMING/{sample}.trimmed.fastq.gz",
    input:
        fastq=lambda wildcards: samples.loc[wildcards.sample].fastq,
    params:
        adapters=" ".join([f"-a {x}" for x in ADAPTERS]),
    log:
        "logs/CUTADAPT_ADAPTER_TRIMMING/{sample}.log",
    conda:
        "../envs/cutadapt.yaml"
    envmodules:
        *config["envmodules"]["cutadapt"],
    benchmark:
        "benchmarks/CUTADAPT_ADAPTER_TRIMMING/{sample}.benchmark.txt"
    message:
        "Cutadapt_Adapter_Trimming: TRIMMING ADAPTERS FOR SAMPLE {input.fastq} WITH CUTADAPT"
    threads: 1
    shell:
        "cutadapt {params.adapters} --minimum-length 31 -o {output.fastq} {input.fastq} &> {log}"


rule FastQC_AfterTrimming:
    output:
        html="results/FASTQC_AFTER_TRIMMING/{sample}.trimmed_fastqc.html",
        zip="results/FASTQC_AFTER_TRIMMING/{sample}.trimmed_fastqc.zip",
    input:
        fastq="results/CUTADAPT_ADAPTER_TRIMMING/{sample}.trimmed.fastq.gz",
    threads: 2
    log:
        "logs/FASTQC_AFTER_TRIMMING/{sample}.log",
    conda:
        "../envs/fastqc.yaml"
    envmodules:
        *config["envmodules"]["fastqc"],
    benchmark:
        "benchmarks/FASTQC_AFTER_TRIMMING/{sample}.benchmark.txt"
    resources:
        mem_mb=10240,
    message:
        "FastQC_AfterTrimming: RUNNING QUALITY CONTROL WITH FASTQC FOR SAMPLE {input.fastq} AFTER TRIMMING ADAPTERS"
    shell:
        "fastqc {input.fastq} --memory {resources.mem_mb} --threads {threads} --nogroup --outdir results/FASTQC_AFTER_TRIMMING &> {log}"


rule MultiQC:
    """Run MultiQC"""
    output:
        html="results/MULTIQC/multiqc_report.html",
    input:
        unpack(multiqc_input),
    log:
        "logs/MULTIQC/MULTIQC.log",
    threads: 1
    conda:
        "../envs/multiqc.yaml"
    params:
        config=os.path.join(WORKFLOW_DIR, "envs", "multiqc_config.yaml"),
    envmodules:
        *config["envmodules"]["multiqc"],
    benchmark:
        "benchmarks/MULTIQC/MULTIQC.benchmark.txt"
    message:
        "MultiQC: COMBINING QUALITY CONTROL METRICS WITH MULTIQC"
    shell:
        'echo {input} | tr " " "\n" > {output.html}.fof;'
        "multiqc -c {params.config} -l {output.html}.fof --verbose --force --outdir results/MULTIQC &> {log}"
