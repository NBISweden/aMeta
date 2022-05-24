checkpoint Create_Sample_TaxID_Directories:
    """Create taxid directory

    For a sample, create taxid for each entry in krakenuniq output
    taxID.pathogens. Downstream rules use the taxid directories as
    input, but it is not known beforehand which these are; they are
    determined by the finds in krakenuniq.

    """
    input:
        pathogens="results/KRAKENUNIQ/{sample}/taxID.pathogens",
    output:
        done="results/AUTHENTICATION/{sample}/.extract_taxids_done",
    log:
        "logs/CREATE_SAMPLE_TAXID_DIRECTORIES/{sample}.log",
    params:
        dir=lambda wildcards: f"results/AUTHENTICATION/{wildcards.sample}",
    shell:
        "mkdir -p {params.dir}; "
        "while read taxid; do mkdir {params.dir}/$taxid; touch {params.dir}/$taxid/.done; done<{input.pathogens};"
        "touch {output.done}"


rule aggregate:
    """aggregate rule: generate all sample/taxid/refid combinations to
    generate targets.

    The reference id depends on rule MaltExtract having been run,
    which therefore needs to be triggered before the aggregate step.
    However, beginning from an empty folder structure, the MaltExtract
    output does not exist which means taxid cannot be converted to
    refid. For this reason, also MaltExtract is a checkpoint, which
    upon completion triggers reevaluation of the DAG.

    """
    input:
        aggregate_maltextract,
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
    """Generate a list of species names for a taxonomic identifier"""
    input:
        dirdone="results/AUTHENTICATION/{sample}/{taxid}/.done",
    output:
        node_list="results/AUTHENTICATION/{sample}/{taxid}/node_list.txt",
    params:
        tax_db=config["krakenuniq_db"],
    log:
        "logs/MAKE_NODE_LIST/{sample}_{taxid}.log",
    shell:
        "awk -v var={wildcards.taxid} '{{ if($1==var) print $0 }}' {params.tax_db}/taxDB | cut -f3 > {output.node_list}"


checkpoint Malt_Extract:
    """Convert rma6 output to misc usable formats.

    Conversion of taxid to refid requires MaltExtract having been run.
    Therefore this rule is a checkpoint that will trigger reevaluation
    of downstream rules. The aggregation is performed by
    aggregate_maltextract.

    """
    input:
        rma6="results/MALT/{sample}.trimmed.rma6",
        node_list="results/AUTHENTICATION/{sample}/{taxid}/node_list.txt",
    output:
        maltextractlog="results/AUTHENTICATION/{sample}/{taxid}/{sample}.trimmed.rma6_MaltExtract_output/log.txt",
        nodeentries="results/AUTHENTICATION/{sample}/{taxid}/{sample}.trimmed.rma6_MaltExtract_output/default/readDist/{sample}.trimmed.rma6_additionalNodeEntries.txt",
    params:
        ncbi_db=config["ncbi_db"],
        extract=format_maltextract_output_directory,
    threads: 4
    log:
        "logs/MALT_EXTRACT/{sample}_{taxid}.log",
    conda:
        "../envs/malt.yaml"
    envmodules:
        *config["envmodules"]["malt"],
    message:
        "Malt_Extract: RUNNING MALT EXTRACT FOR SAMPLE {input.rma6}"
    shell:
        "time MaltExtract -i {input.rma6} -f def_anc -o {params.extract} --reads --threads {threads} --matches --minPI 85.0 --maxReadLength 0 --minComp 0.0 --meganSummary -t {input.node_list} -v 2> {log}"


rule Post_Processing:
    input:
        nodeentries="results/AUTHENTICATION/{sample}/{taxid}/{sample}.trimmed.rma6_MaltExtract_output/default/readDist/{sample}.trimmed.rma6_additionalNodeEntries.txt",
        node_list="results/AUTHENTICATION/{sample}/{taxid}/node_list.txt",
    output:
        analysis="results/AUTHENTICATION/{sample}/{taxid}/{sample}.trimmed.rma6_MaltExtract_output/analysis.RData",
    threads: 4
    params:
        extract=format_maltextract_output_directory,
    log:
        "logs/POST_PROCESSING/{sample}_{taxid}.log",
    conda:
        "../envs/malt.yaml"
    message:
        "Post_Processing: POSTPROCESSING SAMPLES"
    envmodules:
        *config["envmodules"]["malt"],
    shell:
        "postprocessing.AMPS.r -m def_anc -r {params.extract} -t {threads} -n {input.node_list} || echo 'postprocessing failed for {wildcards.sample}_{wildcards.taxid}' > {output.analysis}  2> {log}"


rule Breadth_Of_Coverage:
    input:
        nodeentries="results/AUTHENTICATION/{sample}/{taxid}/{sample}.trimmed.rma6_MaltExtract_output/default/readDist/{sample}.trimmed.rma6_additionalNodeEntries.txt",
        sam="results/MALT/{sample}.trimmed.sam.gz",
    output:
        name_list="results/AUTHENTICATION/{sample}/{taxid}/{refid}/name.list",
        sam="results/AUTHENTICATION/{sample}/{taxid}/{refid}/{taxid}.sam",
        bam="results/AUTHENTICATION/{sample}/{taxid}/{refid}/{taxid}.bam",
        sorted_bam="results/AUTHENTICATION/{sample}/{taxid}/{refid}/{taxid}.sorted.bam",
        breadth_of_coverage="results/AUTHENTICATION/{sample}/{taxid}/{refid}/{taxid}.breadth_of_coverage",
        fasta="results/AUTHENTICATION/{sample}/{taxid}/{refid}/{taxid}.fasta",
    params:
        malt_fasta=config["malt_nt_fasta"],
        ref_id=get_ref_id,
    message:
        "Breadth_Of_Coverage: COMPUTING BREADTH OF COVERAGE, EXTRACTING REFERENCE SEQUENCE FOR VISUALIZING ALIGNMENTS WITH IGV"
    log:
        "logs/BREADTH_OF_COVERAGE/{sample}_{taxid}_{refid}.log",
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
        nodeentries="results/AUTHENTICATION/{sample}/{taxid}/{sample}.trimmed.rma6_MaltExtract_output/default/readDist/{sample}.trimmed.rma6_additionalNodeEntries.txt",
        bam="results/AUTHENTICATION/{sample}/{taxid}/{refid}/{taxid}.sorted.bam",
    output:
        distribution="results/AUTHENTICATION/{sample}/{taxid}/{refid}/{taxid}.read_length.txt",
    message:
        "Read_Length_Distribution: COMPUTING READ LENGTH DISTRIBUTION"
    log:
        "logs/READ_LENGTH_DISTRIBUTION/{sample}_{taxid}_{refid}.log",
    conda:
        "../envs/malt.yaml"
    envmodules:
        *config["envmodules"]["malt"],
    shell:
        "samtools view {input.bam} | awk '{{ print length($10) }}' > {output.distribution}"


rule PMD_scores:
    input:
        nodeentries="results/AUTHENTICATION/{sample}/{taxid}/{sample}.trimmed.rma6_MaltExtract_output/default/readDist/{sample}.trimmed.rma6_additionalNodeEntries.txt",
        bam="results/AUTHENTICATION/{sample}/{taxid}/{refid}/{taxid}.sorted.bam",
    output:
        scores="results/AUTHENTICATION/{sample}/{taxid}/{refid}/{taxid}.PMDscores.txt",
    message:
        "PMD_scores: COMPUTING PMD SCORES"
    log:
        "logs/PMD_SCORES/{sample}_{taxid}_{refid}.log",
    conda:
        "../envs/malt.yaml"
    envmodules:
        *config["envmodules"]["malt"],
    shell:
        "samtools view -h {input.bam} | pmdtools --number 100000 --printDS > {output.scores}"


rule Authentication_Plots:
    input:
        dir="results/AUTHENTICATION/{sample}/{taxid}",
        nodeentries="results/AUTHENTICATION/{sample}/{taxid}/{sample}.trimmed.rma6_MaltExtract_output/default/readDist/{sample}.trimmed.rma6_additionalNodeEntries.txt",
        node_list="results/AUTHENTICATION/{sample}/{taxid}/node_list.txt",
        distribution="results/AUTHENTICATION/{sample}/{taxid}/{refid}/{taxid}.read_length.txt",
        scores="results/AUTHENTICATION/{sample}/{taxid}/{refid}/{taxid}.PMDscores.txt",
        breadth_of_coverage="results/AUTHENTICATION/{sample}/{taxid}/{refid}/{taxid}.breadth_of_coverage",
    output:
        plot="results/AUTHENTICATION/{sample}/{taxid}/{refid}/authentic_Sample_{sample}.trimmed.rma6_TaxID_{taxid}.pdf",
    params:
        exe=WORKFLOW_DIR / "scripts/authentic.R",
    message:
        "Authentication_Plots: MAKING AUTHENTICATION AND VALIDATION PLOTS"
    log:
        "logs/AUTHENTICATION_PLOTS/{sample}_{taxid}_{refid}.log",
    conda:
        "../envs/malt.yaml"
    envmodules:
        *config["envmodules"]["malt"],
    shell:
        "Rscript {params.exe} {wildcards.taxid} {wildcards.sample}.trimmed.rma6 {input.dir}/{wildcards.refid}"


rule Deamination:
    input:
        bam="results/AUTHENTICATION/{sample}/{taxid}/{refid}/{taxid}.sorted.bam",
    output:
        tmp="results/AUTHENTICATION/{sample}/{taxid}/{refid}/PMD_temp.txt",
        pmd="results/AUTHENTICATION/{sample}/{taxid}/{refid}/PMD_plot.frag.pdf",
    message:
        "Deamination: INFERRING DEAMINATION PATTERN FROM CPG SITES"
    log:
        "logs/DEAMINATION/{sample}_{taxid}_{refid}.log",
    conda:
        "../envs/malt.yaml"
    envmodules:
        *config["envmodules"]["malt"],
    shell:
        "samtools view {input.bam} | pmdtools --platypus > {output.tmp}; "
        "cd results/AUTHENTICATION/{wildcards.sample}/{wildcards.taxid}/{wildcards.refid}; "
        "R CMD BATCH $(which plotPMD); "
