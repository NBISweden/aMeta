rule Build_Malt_DB:
    output:
        seqid2taxid_project="results/MALT_DB/seqid2taxid.project.map",
        seqids_project="results/MALT_DB/seqids.project",
        project_headers="results/MALT_DB/project.headers",
        project_fasta="results/MALT_DB/library.project.fna",
        db=directory("results/MALT_DB/maltDB.dat"),
    input:
        unique_taxids="results/KRAKENUNIQ_ABUNDANCE_MATRIX/unique_species_taxid_list.txt",
    params:
        seqid2taxid=config["malt_seqid2taxid_db"],
        nt_fasta=config["malt_nt_fasta"],
        accession2taxid=config["malt_accession2taxid"],
    threads: 20
    log:
        "logs/BUILD_MALT_DB/BUILD_MALT_DB.log",
    conda:
        "../envs/malt.yaml"
    envmodules:
        *config["envmodules"]["Build_Malt_DB"],
    benchmark:
        "benchmarks/BUILD_MALT_DB/BUILD_MALT_DB.benchmark.txt"
    message:
        "BUILDING MALT DATABASE USING SPECIES DETECTED BY KRAKENUNIQ"
    shell:
        "grep -wFf {input.unique_taxids} {params.seqid2taxid} > {output.seqid2taxid_project}; "
        "cut -f1 {output.seqid2taxid_project} > {output.seqids_project}; "
        "grep -Ff {output.seqids_project} {params.nt_fasta} | sed 's/>//g' > {output.project_headers}; "
        "seqtk subseq {params.nt_fasta} {output.project_headers} > {output.project_fasta} 2> {log}; "
        "unset DISPLAY; "
        "malt-build -i {output.project_fasta} -a2taxonomy {params.accession2taxid} -s DNA -t {threads} -d {output.db} &>> {log}"


rule Malt:
    output:
        rma6="results/MALT/{sample}.trimmed.rma6",
        sam="results/MALT/{sample}.trimmed.sam.gz",
    input:
        fastq="results/CUTADAPT_ADAPTER_TRIMMING/{sample}.trimmed.fastq.gz",
        db="results/MALT_DB/maltDB.dat",
    params:
        gunzipped_sam="results/MALT/{sample}.trimmed.sam",
    threads: 20
    log:
        "logs/MALT/{sample}.log",
    conda:
        "../envs/malt.yaml"
    envmodules:
        *config["envmodules"]["Malt"],
    benchmark:
        "benchmarks/MALT/{sample}.benchmark.txt"
    message:
        "RUNNING MALT ALIGNMENTS FOR SAMPLE {input.fastq}"
    shell:
        "unset DISPLAY; malt-run -at SemiGlobal -m BlastN -i {input.fastq} -o {output.rma6} -a {params.gunzipped_sam} -t {threads} -d {input.db} &> {log}"


rule Malt_QuantifyAbundance:
    output:
        out_dir=directory("results/MALT_QUANTIFY_ABUNDANCE/{sample}"),
        counts="results/MALT_QUANTIFY_ABUNDANCE/{sample}/sam_counts.txt",
    input:
        sam="results/MALT/{sample}.trimmed.sam.gz",
    params:
        unique_taxids="results/KRAKENUNIQ_ABUNDANCE_MATRIX/unique_species_taxid_list.txt",
    log:
        "logs/MALT_QUANTIFY_ABUNDANCE/{sample}.log",
    benchmark:
        "benchmarks/MALT_QUANTIFY_ABUNDANCE/{sample}.benchmark.txt"
    message:
        "QUANTIFYING MICROBIAL ABUNDANCE USING MALT ALIGNMENTS FOR SAMPLE {input.sam}"
    shell:
        """n_species=0; touch {output.counts}; for i in $(cat {params.unique_taxids}); do zgrep \"|tax|$i\" {input.sam} | grep -v '@' | awk '!a[$1]++' | wc -l >> {output.counts} || true; n_species=$((n_species+1)); echo Finished $n_species species; done &> {log}"""


rule Malt_AbundanceMatrix:
    output:
        out_dir=directory("results/MALT_ABUNDANCE_MATRIX"),
        abundance_matrix="results/MALT_ABUNDANCE_MATRIX/malt_abundance_matrix.txt",
    input:
        sam_counts=expand(
            "results/MALT_QUANTIFY_ABUNDANCE/{sample}/sam_counts.txt", sample=SAMPLES
        ),
    log:
        "logs/MALT_ABUNDANCE_MATRIX/MALT_ABUNDANCE_MATRIX.log",
    params:
        exe=WORKFLOW_DIR / "scripts/malt_abundance_matrix.R",
    conda:
        "../envs/r.yaml"
    envmodules:
        *config["envmodules"]["Malt_AbundanceMatrix"],
    benchmark:
        "benchmarks/MALT_ABUNDANCE_MATRIX/MALT_ABUNDANCE_MATRIX.benchmark.txt"
    message:
        "COMPUTING MALT MICROBIAL ABUNDANCE MATRIX"
    shell:
        "Rscript {params.exe} results/MALT_QUANTIFY_ABUNDANCE {output.out_dir} &> {log}"


rule Authentication:
    output:
        out_dir=directory("results/AUTHENTICATION/{sample}"),
    input:
        rma6="results/MALT/{sample}.trimmed.rma6",
        sam="results/MALT/{sample}.trimmed.sam.gz",
        pathogen_tax_id="results/KRAKENUNIQ/{sample}/taxID.pathogens",
    params:
        krakenuniq_db=config["krakenuniq_db"],
        ncbi_db=config["ncbi_db"],
        malt_fasta=config["malt_nt_fasta"],
        exe=WORKFLOW_DIR / "scripts/authentic.sh",
        scripts_dir=WORKFLOW_DIR / "scripts"
    log:
        "logs/AUTHENTICATION/{sample}.AUTHENTICATION.log",
    conda:
        "../envs/malt.yaml"
    benchmark:
        "benchmarks/AUTHENTICATION/{sample}.AUTHENTICATION.benchmark.txt"
    message:
        "PERFORMING AUTHENTICATION ANALYSIS ON MALT ALIGNMENTS FOR SAMPLE {input.rma6}"
    shell:
        "mkdir -p results/AUTHENTICATION || true &> {log}; "
        "mkdir {output.out_dir} || true &> {log}; "
        "for i in $(cat {input.pathogen_tax_id}); do echo Authenticating taxon $i; {params.exe} $i results/MALT {input.rma6} {input.sam} {output.out_dir}/$i {params.scripts_dir} {params.krakenuniq_db} {params.ncbi_db} {params.malt_fasta}; done &> {log}"
