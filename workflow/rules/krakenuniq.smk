rule KrakenUniq:
    """Run KrakenUniq on trimmed fastq data"""
    output:
        report = "results/KRAKENUNIQ/{sample}/krakenuniq.output",
        seqs = "results/KRAKENUNIQ/{sample}/sequences.krakenuniq"
    input:
        fastq = "results/CUTADAPT_ADAPTER_TRIMMING/{sample}.trimmed.fastq.gz"
    params:
        #DB = "/proj/snic2018-8-150/uppstore2018095/private/NBIS_Demo/DBDIR_KrakenUniq_MicrobialNT_Plus_CompleteGenomes"
        DB = config["krakenuniq_db"]
    threads: 10
    log:
        "logs/KRAKENUNIQ/{sample}.log"
    conda:
        "../envs/krakenuniq.yaml"
    envmodules:
        *config["envmodules"]["KrakenUniq"],
    benchmark:
        "benchmarks/KRAKENUNIQ/{sample}.benchmark.txt"
    message:
        "PERFORMING TAXONOMIC CLASSIFICATION OF SAMPLE {input.fastq} WITH KRAKENUNIQ"
    shell:
        "krakenuniq --db {params.DB} --fastq-input {input.fastq} --threads {threads} --output {output.seqs} --report-file {output.report} --gzip-compressed --only-classified-out &> {log}"


rule Filter_KrakenUniq_Output:
    output:
        filtered = "results/KRAKENUNIQ/{sample}/krakenuniq.output.filtered",
        pathogens = "results/KRAKENUNIQ/{sample}/krakenuniq.output.pathogens",
        pathogen_tax_id = "results/KRAKENUNIQ/{sample}/taxID.pathogens"
    input:
        "results/KRAKENUNIQ/{sample}/krakenuniq.output"
    log:
        "logs/FILTER_KRAKENUNIQ_OUTPUT/{sample}.log"
    benchmark:
        "benchmarks/FILTER_KRAKENUNIQ_OUTPUT/{sample}.benchmark.txt"
    message:
        "APPLYING DEPTH AND BREADTH OF COVERAGE FILTERS TO KRAKENUNIQ OUTPUT FOR SAMPLE {input}"
    shell:
        """Rscript scripts/filter_krakenuniq.R {input} &> {log}"""
        """cut -f7 {output.pathogens} | tail -n +2 > {output.pathogen_tax_id}"""


rule KrakenUniq2Krona:
    output:
        tax_ids = "results/KRAKENUNIQ/{sample}/krakenuniq.output_taxIDs_kmers1000.txt",
        seqs = "results/KRAKENUNIQ/{sample}/sequences.krakenuniq_kmers1000.txt",
        krona = "results/KRAKENUNIQ/{sample}/sequences.krakenuniq_kmers1000.krona",
        html = "results/KRAKENUNIQ/{sample}/taxonomy.krona.html"
    input:
        report = "results/KRAKENUNIQ/{sample}/krakenuniq.output",
        seqs = "results/KRAKENUNIQ/{sample}/sequences.krakenuniq"
    log:
        "logs/KRAKENUNIQ2KRONA/{sample}.log"
    conda:
        "../envs/krona.yaml"
    envmodules:
        *config["envmodules"]["KrakenUniq2Krona"],
    benchmark:
        "benchmarks/KRAKENUNIQ2KRONA/{sample}.benchmark.txt"
    message:
        "VISUALIZING KRAKENUNIQ RESULTS WITH KRONA FOR SAMPLE {input.report}"
    shell:
        "Rscript scripts/krakenuniq2krona.R {input.report} {input.seqs} &> {log}; "
        "cat {output.seqs} | cut -f 2,3 > {output.krona}; "
        "ktImportTaxonomy {output.krona} -o {output.html} &>> {log}"


rule KrakenUniq_AbundanceMatrix:
    output:
        out_dir = directory("results/KRAKENUNIQ_ABUNDANCE_MATRIX"),
        unique_species = "results/KRAKENUNIQ_ABUNDANCE_MATRIX/unique_species_taxid_list.txt",
        unique_species_names = "results/KRAKENUNIQ_ABUNDANCE_MATRIX/unique_species_names_list.txt",
        abundance_matrix = "results/KRAKENUNIQ_ABUNDANCE_MATRIX/krakenuniq_abundance_matrix.txt"
    input:
        expand("results/KRAKENUNIQ/{sample}/krakenuniq.output", sample = SAMPLES)
    log:
        "logs/KRAKENUNIQ_ABUNDANCE_MATRIX/KRAKENUNIQ_ABUNDANCE_MATRIX.log"
    benchmark:
        "benchmarks/KRAKENUNIQ_ABUNDANCE_MATRIX/KRAKENUNIQ_ABUNDANCE_MATRIX.benchmark.txt"
    message:
        "COMPUTING KRAKENUNIQ MICROBIAL ABUNDANCE MATRIX"
    shell:
        "Rscript scripts/krakenuniq_abundance_matrix.R results/KRAKENUNIQ {output.out_dir} &> {log}"
