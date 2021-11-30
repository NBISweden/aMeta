rule FastQC_BeforeTrimming:
    """Run fastq before trimming"""
    output:
        html="results/FASTQC_BEFORE_TRIMMING/{sample}_fastqc.html",
    input:
        fastq="data/{sample}.fastq.gz",
    params:
        outdir="results/FASTQC_BEFORE_TRIMMING",
    conda:
        "../envs/fastqc.yaml"
    envmodules:
        *config["envmodules"]["FastQC_BeforeTrimming"],
    log:
        "logs/FASTQC_BEFORE_TRIMMING/{sample}.log",
    threads: 1
    shell:
        "fastqc {input.fastq} --threads {threads} --nogroup --outdir {params.outdir} &> {log}"


rule MultiQC_BeforeTrimming:
        input:
                expand("results/FASTQC_BEFORE_TRIMMING/{sample}_fastqc.html", sample = SAMPLES)
        params:
                out_dir = "results/MULTIQC_BEFORE_TRIMMING"
        log:
                "logs/MULTIQC_BEFORE_TRIMMING/MULTIQC_BEFORE_TRIMMING.log"
        benchmark:
                "benchmarks/MULTIQC_BEFORE_TRIMMING/MULTIQC_BEFORE_TRIMMING.benchmark.txt"
        message:
                "COMBINING QUALITY CONTROL METRICS WITH MULTIQC BEFORE TRIMMING ADAPTERS"
        output:
                html = "results/MULTIQC_BEFORE_TRIMMING/multiqc_report.html"
        run:
                shell("(multiqc results/FASTQC_BEFORE_TRIMMING --verbose --force --outdir {params.out_dir}) &> {log}")
