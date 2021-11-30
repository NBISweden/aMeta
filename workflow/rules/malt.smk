rule Build_Malt_DB:
    output:
        seqid2taxid_project = "results/MALT_DB/seqid2taxid.project.map",
        seqids_project = "results/MALT_DB/seqids.project",
        project_headers = "results/MALT_DB/project.headers",
        project_fasta = "results/MALT_DB/library.project.fna",
        db = directory("results/MALT_DB/maltDB.dat")
    input:
        unique_taxids = "results/KRAKENUNIQ_ABUNDANCE_MATRIX/unique_species_taxid_list.txt"
    params:
        seqid2taxid = "/proj/snic2018-8-150/uppstore2018095/private/NBIS_Demo/DBDIR_KrakenUniq_Full_NT/seqid2taxid.map.orig",
        nt_fasta = "/proj/snic2018-8-150/uppstore2018095/private/NBIS_Demo/DBDIR_KrakenUniq_Full_NT/library/nt/library.fna",
        accession2taxid = "/proj/snic2018-8-150/uppstore2018095/private/NBIS_Demo/DBDIR_KrakenUniq_Full_NT/taxonomy/nucl_gb.accession2taxid"
    threads: 20
    log:
        "logs/BUILD_MALT_DB/BUILD_MALT_DB.log"
    benchmark:
        "benchmarks/BUILD_MALT_DB/BUILD_MALT_DB.benchmark.txt"
    message:
        "BUILDING MALT DATABASE USING SPECIES DETECTED BY KRAKENUNIQ"
    run:
        shell("grep -wFf {input.unique_taxids} {params.seqid2taxid} > {output.seqid2taxid_project}")
        shell("cut -f1 {output.seqid2taxid_project} > {output.seqids_project}")
        shell("grep -Ff {output.seqids_project} {params.nt_fasta} | sed 's/>//g' > {output.project_headers}")
        shell("(seqtk subseq {params.nt_fasta} {output.project_headers} > {output.project_fasta}) 2> {log}")
        shell("unset DISPLAY")
        shell("(" + MALT_HOME + "/malt-build -i {output.project_fasta} -a2taxonomy {params.accession2taxid} -s DNA -t {threads} -d {output.db}) &> {log}")

rule Malt:
        input:
                fastq = "results/CUTADAPT_ADAPTER_TRIMMING/{sample}.trimmed.fastq.gz",
                db = "results/MALT_DB/maltDB.dat"
        params:
                gunzipped_sam = "results/MALT/{sample}.trimmed.sam"
        threads: 20
        log:
                "logs/MALT/{sample}.log"
        benchmark:
                "benchmarks/MALT/{sample}.benchmark.txt"
        message:
                "RUNNING MALT ALIGNMENTS FOR SAMPLE {input.fastq}"
        output:
                rma6 = "results/MALT/{sample}.trimmed.rma6",
                sam = "results/MALT/{sample}.trimmed.sam.gz"
        run:
                shell("unset DISPLAY")
                shell("(" + MALT_HOME + "/malt-run -at SemiGlobal -m BlastN -i {input.fastq} -o {output.rma6} -a {params.gunzipped_sam} -t {threads} -d {input.db}) &> {log}")

rule Malt_QuantifyAbundance:
        input:
                sam = "results/MALT/{sample}.trimmed.sam.gz"
        params:
                unique_taxids = "results/KRAKENUNIQ_ABUNDANCE_MATRIX/unique_species_taxid_list.txt"
        log:
                "logs/MALT_QUANTIFY_ABUNDANCE/{sample}.log"
        benchmark:
                "benchmarks/MALT_QUANTIFY_ABUNDANCE/{sample}.benchmark.txt"
        message:
                "QUANTIFYING MICROBIAL ABUNDANCE USING MALT ALIGNMENTS FOR SAMPLE {input.sam}"
        output:
                out_dir = directory("results/MALT_QUANTIFY_ABUNDANCE/{sample}"),
                counts = "results/MALT_QUANTIFY_ABUNDANCE/{sample}/sam_counts.txt"
        run:
                shell("(n_species=0; touch {output.counts}; for i in $(cat {params.unique_taxids}); do zgrep \"|tax|$i|\" {input.sam} | grep -v '@' | awk '!a[$1]++' | wc -l >> {output.counts} || true; n_species=$((n_species+1)); echo Finished $n_species species; done) &> {log}")

rule Malt_AbundanceMatrix:
        input:
                sam_counts = expand("results/MALT_QUANTIFY_ABUNDANCE/{sample}/sam_counts.txt", sample=SAMPLES)
        log:
                "logs/MALT_ABUNDANCE_MATRIX/MALT_ABUNDANCE_MATRIX.log"
        benchmark:
                "benchmarks/MALT_ABUNDANCE_MATRIX/MALT_ABUNDANCE_MATRIX.benchmark.txt"
        message:
                "COMPUTING MALT MICROBIAL ABUNDANCE MATRIX"
        output:
                out_dir = directory("results/MALT_ABUNDANCE_MATRIX"),
                abundance_matrix = "results/MALT_ABUNDANCE_MATRIX/malt_abundance_matrix.txt"
        run:
                shell("(Rscript scripts/malt_abundance_matrix.R results/MALT_QUANTIFY_ABUNDANCE {output.out_dir}) &> {log}")
