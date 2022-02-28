# the checkpoint that shall trigger re-evaluation of the DAG
checkpoint Extract_TaxIDs:
    input:
        pathogens="results/KRAKENUNIQ/{sample}/taxID.pathogens",
    output:
        dir=directory("results/AUTHENTICATION/{sample}"),
    log:
        "logs/EXTRACT_TAXIDS/{sample}.log",
    shell:
        "mkdir -p {output.dir}; "
        "while read taxid; do mkdir {output.dir}/$taxid; done<{input.pathogens}"


rule aggregate:
    input:
        aggregate_PMD,
        aggregate_plots,
        aggregate_post,
    output:
        "results/AUTHENTICATION/.{sample}_done",
    log:
        "logs/AGGREGATE/{sample}.log",
    shell:
        "touch {output}; "


rule Make_Node_List:
    input:
        dir="results/AUTHENTICATION/{sample}/{taxid}/",
    output:
        node_list="results/AUTHENTICATION/{sample}/{taxid}/node_list.txt",
    params:
        tax_db=config["krakenuniq_db"],
    log:
        "logs/MAKE_NODE_LIST/{sample}_{taxid}.log",
    shell:
        "awk -v var={wildcards.taxid} '{{ if($1==var) print $0 }}' {params.tax_db}/taxDB | cut -f3 > {output.node_list}"


rule Malt_Extract:
    input:
        rma6="results/MALT/{sample}.trimmed.rma6",
        node_list="results/AUTHENTICATION/{sample}/{taxid}/node_list.txt",
    output:
        extract=directory(
            "results/AUTHENTICATION/{sample}/{taxid}/{sample}.trimmed.rma6_MaltExtract_output"
        ),
    params:
        ncbi_db=config["ncbi_db"],
    threads: 4
    log:
        "logs/MALT_EXTRACT/{sample}_{taxid}.log",
    conda:
        "../envs/malt.yaml"
    envmodules:
        *config["envmodules"]["malt"],
    message:
        "RUNNING MALT EXTRACT FOR SAMPLE {input.rma6}"
    shell:
        "time MaltExtract -i {input.rma6} -f def_anc -o {output.extract} --reads --threads {threads} --matches --minPI 85.0 --maxReadLength 0 --minComp 0.0 --meganSummary -r {params.ncbi_db} -t {input.node_list} -v"


rule Post_Processing:
    input:
        extract="results/AUTHENTICATION/{sample}/{taxid}/{sample}.trimmed.rma6_MaltExtract_output",
        node_list="results/AUTHENTICATION/{sample}/{taxid}/node_list.txt",
    output:
        analysis="results/AUTHENTICATION/{sample}/{taxid}/{sample}.trimmed.rma6_MaltExtract_output/analysis.RData",
    threads: 4
    log:
        "logs/POST_PROCESSING/{sample}_{taxid}.log",
    conda:
        "../envs/malt.yaml"
    envmodules:
        *config["envmodules"]["malt"],
    shell:
        "postprocessing.AMPS.r -m def_anc -r {input.extract} -t {threads} -n {input.node_list}"


rule Breadth_Of_Coverage:
    input:
        extract="results/AUTHENTICATION/{sample}/{taxid}/{sample}.trimmed.rma6_MaltExtract_output",
        sam="results/MALT/{sample}.trimmed.sam.gz",
    output:
        name_list="results/AUTHENTICATION/{sample}/{taxid}/name.list",
        sam="results/AUTHENTICATION/{sample}/{taxid}/{taxid}.sam",
        bam="results/AUTHENTICATION/{sample}/{taxid}/{taxid}.bam",
        sorted_bam="results/AUTHENTICATION/{sample}/{taxid}/{taxid}.sorted.bam",
        breadth_of_coverage="results/AUTHENTICATION/{sample}/{taxid}/{taxid}.breadth_of_coverage",
        fasta="results/AUTHENTICATION/{sample}/{taxid}/{taxid}.fasta",
    params:
        malt_fasta=config["malt_nt_fasta"],
        ref_id=get_ref_id,
    message:
        "COMPUTING BREADTH OF COVERAGE, EXTRACTING REFERENCE SEQUENCE FOR VISUALIZING ALIGNMENTS WITH IGV"
    log:
        "logs/BREADTH_OF_COVERAGE/{sample}_{taxid}.log",
    conda:
        "../envs/malt.yaml"
    envmodules:
        *config["envmodules"]["malt"],
    shell:
        "echo {params.ref_id} > {output.name_list}; "
        "zgrep {params.ref_id} {input.sam} | uniq > {output.sam}; "
        "samtools view -bS {output.sam} > {output.bam}; "
        "samtools sort {output.bam} > {output.sorted_bam}; "
        "samtools index {output.sorted_bam}; "
        "samtools depth -a {output.sorted_bam} > {output.breadth_of_coverage}; "
        "seqtk subseq {params.malt_fasta} {output.name_list} > {output.fasta}"


rule Read_Length_Distribution:
    input:
        bam="results/AUTHENTICATION/{sample}/{taxid}/{taxid}.sorted.bam",
    output:
        distribution="results/AUTHENTICATION/{sample}/{taxid}/{taxid}.read_length.txt",
    message:
        "COMPUTING READ LENGTH DISTRIBUTION"
    log:
        "logs/READ_LENGTH_DISTRIBUTION/{sample}_{taxid}.log",
    conda:
        "../envs/malt.yaml"
    envmodules:
        *config["envmodules"]["malt"],
    shell:
        "samtools view {input.bam} | awk '{{ print length($10) }}' > {output.distribution}"


rule PMD_scores:
    input:
        bam="results/AUTHENTICATION/{sample}/{taxid}/{taxid}.sorted.bam",
    output:
        scores="results/AUTHENTICATION/{sample}/{taxid}/{taxid}.PMDscores.txt",
    message:
        "COMPUTING PMD SCORES"
    log:
        "logs/PMD_SCORES/{sample}_{taxid}.log",
    conda:
        "../envs/malt.yaml"
    envmodules:
        *config["envmodules"]["malt"],
    shell:
        "samtools view -h {input.bam} | pmdtools --printDS > {output.scores}"


rule Authentication_Plots:
    input:
        dir="results/AUTHENTICATION/{sample}/{taxid}",
        node_list="results/AUTHENTICATION/{sample}/{taxid}/node_list.txt",
        distribution="results/AUTHENTICATION/{sample}/{taxid}/{taxid}.read_length.txt",
        scores="results/AUTHENTICATION/{sample}/{taxid}/{taxid}.PMDscores.txt",
        breadth_of_coverage="results/AUTHENTICATION/{sample}/{taxid}/{taxid}.breadth_of_coverage",
    output:
        plot="results/AUTHENTICATION/{sample}/{taxid}/authentic_Sample_{sample}.trimmed.rma6_TaxID_{taxid}.pdf",
    params:
        exe=WORKFLOW_DIR / "scripts/authentic.R",
    message:
        "MAKING AUTHENTICATION AND VALIDATION PLOTS"
    log:
        "logs/AUTHENTICATION_PLOTS/{sample}_{taxid}.log",
    conda:
        "../envs/malt.yaml"
    envmodules:
        *config["envmodules"]["malt"],
    shell:
        "Rscript {params.exe} {wildcards.taxid} {wildcards.sample}.trimmed.rma6 {input.dir}"


rule Deamination:
    input:
        bam="results/AUTHENTICATION/{sample}/{taxid}/{taxid}.sorted.bam",
    output:
        tmp="results/AUTHENTICATION/{sample}/{taxid}/PMD_temp.txt",
        pmd="results/AUTHENTICATION/{sample}/{taxid}/PMD_plot.frag.pdf",
    message:
        "INFERRING DEAMINATION PATTERN FROM CPG SITES"
    log:
        "logs/DEAMINATION/{sample}_{taxid}.log",
    conda:
        "../envs/malt.yaml"
    envmodules:
        *config["envmodules"]["malt"],
    shell:
        "samtools view {input.bam} | pmdtools --platypus > {output.tmp}; "
        "cd results/AUTHENTICATION/{wildcards.sample}/{wildcards.taxid}; "
        "R CMD BATCH $(which plotPMD); "
