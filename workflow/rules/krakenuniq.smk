rule KrakenUniq:
    """Run KrakenUniq on trimmed fastq data"""
    output:
        report="results/KRAKENUNIQ/{sample}/krakenuniq.output",
        seqs="results/KRAKENUNIQ/{sample}/sequences.krakenuniq",
    input:
        fastq="results/CUTADAPT_ADAPTER_TRIMMING/{sample}.trimmed.fastq.gz",
    params:
        DB=config["krakenuniq_db"],
    threads: 10
    log:
        "logs/KRAKENUNIQ/{sample}.log",
    conda:
        "../envs/krakenuniq.yaml"
    envmodules:
        *config["envmodules"]["krakenuniq"],
    benchmark:
        "benchmarks/KRAKENUNIQ/{sample}.benchmark.txt"
    message:
        "KrakenUniq: PERFORMING TAXONOMIC CLASSIFICATION OF SAMPLE {input.fastq} WITH KRAKENUNIQ"
    shell:
        "krakenuniq --preload --db {params.DB} --fastq-input {input.fastq} --threads {threads} --output {output.seqs} --report-file {output.report} --gzip-compressed --only-classified-out &> {log}"


rule Filter_KrakenUniq_Output:
    output:
        filtered="results/KRAKENUNIQ/{sample}/krakenuniq.output.filtered",
        pathogens="results/KRAKENUNIQ/{sample}/krakenuniq.output.pathogens",
        pathogen_tax_id="results/KRAKENUNIQ/{sample}/taxID.pathogens",
    input:
        krakenuniq="results/KRAKENUNIQ/{sample}/krakenuniq.output",
        pathogenomesFound=config["pathogenomesFound"],
    log:
        "logs/FILTER_KRAKENUNIQ_OUTPUT/{sample}.log",
    params:
        exe=WORKFLOW_DIR / "scripts/filter_krakenuniq.py",
        n_unique_kmers=config["n_unique_kmers"],
        n_tax_reads=config["n_tax_reads"],
    conda:
        "../envs/krakenuniq.yaml"
    envmodules:
        *config["envmodules"]["krakenuniq"],
    benchmark:
        "benchmarks/FILTER_KRAKENUNIQ_OUTPUT/{sample}.benchmark.txt"
    message:
        "Filter_KrakenUniq_Output: APPLYING DEPTH AND BREADTH OF COVERAGE FILTERS TO KRAKENUNIQ OUTPUT FOR SAMPLE {input}"
    shell:
        """{params.exe} {input.krakenuniq} {params.n_unique_kmers} {params.n_tax_reads} {input.pathogenomesFound} &> {log}; """
        """cut -f7 {output.pathogens} | tail -n +2 > {output.pathogen_tax_id}"""


rule KrakenUniq2Krona:
    output:
        tax_ids="results/KRAKENUNIQ/{sample}/krakenuniq.output.filtered_taxIDs_kmers1000.txt",
        seqs="results/KRAKENUNIQ/{sample}/sequences.krakenuniq_kmers1000.txt",
        krona="results/KRAKENUNIQ/{sample}/sequences.krakenuniq_kmers1000.krona",
        html="results/KRAKENUNIQ/{sample}/taxonomy.krona.html",
    input:
        report="results/KRAKENUNIQ/{sample}/krakenuniq.output.filtered",
        seqs="results/KRAKENUNIQ/{sample}/sequences.krakenuniq",
    log:
        "logs/KRAKENUNIQ2KRONA/{sample}.log",
    conda:
        "../envs/krona.yaml"
    envmodules:
        *config["envmodules"]["krona"],
    params:
        exe=WORKFLOW_DIR / "scripts/krakenuniq2krona.py",
        DB=f"--tax {config['krona_db']}" if "krona_db" in config else "",
    benchmark:
        "benchmarks/KRAKENUNIQ2KRONA/{sample}.benchmark.txt"
    message:
        "KrakenUniq2Krona: VISUALIZING KRAKENUNIQ RESULTS WITH KRONA FOR SAMPLE {input.report}"
    shell:
        "{params.exe} {input.report} {input.seqs} &> {log}; "
        "cat {output.seqs} | cut -f 2,3 > {output.krona}; "
        "ktImportTaxonomy {output.krona} -o {output.html} {params.DB} &>> {log}"


rule KrakenUniq_AbundanceMatrix:
    output:
        out_dir=directory("results/KRAKENUNIQ_ABUNDANCE_MATRIX"),
        unique_species="results/KRAKENUNIQ_ABUNDANCE_MATRIX/unique_species_taxid_list.txt",
        unique_species_names="results/KRAKENUNIQ_ABUNDANCE_MATRIX/unique_species_names_list.txt",
        abundance_matrix="results/KRAKENUNIQ_ABUNDANCE_MATRIX/krakenuniq_abundance_matrix.txt",
    input:
        expand("results/KRAKENUNIQ/{sample}/krakenuniq.output.filtered", sample=SAMPLES),
    log:
        "logs/KRAKENUNIQ_ABUNDANCE_MATRIX/KRAKENUNIQ_ABUNDANCE_MATRIX.log",
    params:
        exe=WORKFLOW_DIR / "scripts/krakenuniq_abundance_matrix.R",
        n_unique_kmers=config["n_unique_kmers"],
        n_tax_reads=config["n_tax_reads"],
    conda:
        "../envs/r.yaml"
    envmodules:
        *config["envmodules"]["r"],
    benchmark:
        "benchmarks/KRAKENUNIQ_ABUNDANCE_MATRIX/KRAKENUNIQ_ABUNDANCE_MATRIX.benchmark.txt"
    message:
        "KrakenUniq_AbundanceMatrix: COMPUTING KRAKENUNIQ MICROBIAL ABUNDANCE MATRIX"
    shell:
        "Rscript {params.exe} results/KRAKENUNIQ {output.out_dir} {params.n_unique_kmers} {params.n_tax_reads} &> {log}"
