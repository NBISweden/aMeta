
# the checkpoint that shall trigger re-evaluation of the DAG
checkpoint Extract_TaxIDs:
    input:
        pathogens="results/KRAKENUNIQ/{sample}/taxID.pathogens"
    output:
        dir=directory("results/AUTHENTICATION/{sample}")
    shell:
        "mkdir -p results/AUTHENTICATION/{wildcards.sample}; "
        "while read taxid; do mkdir results/AUTHENTICATION/{wildcards.sample}/$taxid; done<{input.pathogens}"

def aggregate_PMD(wildcards):
    checkpoint_output = checkpoints.Extract_TaxIDs.get(sample=wildcards.sample).output[0]
    return expand("results/AUTHENTICATION/{sample}/{taxid}/PMD_temp.txt", sample=wildcards.sample, taxid=glob_wildcards(os.path.join(checkpoint_output, "{taxid,[0-9]+}")).taxid)

def aggregate_plots(wildcards):
    checkpoint_output = checkpoints.Extract_TaxIDs.get(sample=wildcards.sample).output[0]
    return expand("results/AUTHENTICATION/{sample}/{taxid}/authentic_Sample_{sample}.trimmed.rma6_TaxID_{taxid}.pdf", sample=wildcards.sample, taxid=glob_wildcards(os.path.join(checkpoint_output, "{taxid,[0-9]+}")).taxid)

def aggregate_post(wildcards):
    checkpoint_output = checkpoints.Extract_TaxIDs.get(sample=wildcards.sample).output[0]
    return expand("results/AUTHENTICATION/{sample}/{taxid}/{sample}.trimmed.rma6_MaltExtract_output/analysis.RData", sample=wildcards.sample, taxid=glob_wildcards(os.path.join(checkpoint_output, "{taxid,[0-9]+}")).taxid)

rule aggregate:
    input:
        aggregate_PMD,
        aggregate_plots,
        aggregate_post
    output:
        "results/AUTHENTICATION/{sample}_status/done"
    shell:
        "mkdir -p results/AUTHENTICATION/{wildcards.sample}_status/; "
        "touch {output}"

#awk -v var="$TAXID" '{if($1==var)print$0}' $TAXDB_DIR/taxDB | cut -f3 > $OUT_DIR/node_list.txt
rule Make_Node_List:
    input:
        dir="results/AUTHENTICATION/{sample}/{taxid}/"
    output:
        node_list="results/AUTHENTICATION/{sample}/{taxid,[0-9]+}/node_list.txt"
    params:
        tax_db=config["krakenuniq_db"]
    shell:
        "TAXID=$(basename {input.dir});"
        "awk -v var=\"$TAXID\" '{{ if($1==var) print $0 }}' {params.tax_db}/taxDB | cut -f3 > {output.node_list}"
        
        
# time MaltExtract -i $IN_DIR/$RMA6 -f def_anc -o $OUT_DIR/${RMA6}_MaltExtract_output --reads --threads $THREADS --matches --minPI 85.0 --maxReadLength 0 --minComp 0.0 --meganSummary -r $NCBI_DB -t $OUT_DIR/node_list.txt -v
rule Malt_Extract:
    input:
        rma6="results/MALT/{sample}.trimmed.rma6",
        node_list="results/AUTHENTICATION/{sample}/{taxid}/node_list.txt"
    output:
        extract=directory("results/AUTHENTICATION/{sample}/{taxid,[0-9]+}/{sample}.trimmed.rma6_MaltExtract_output")
    params:
        ncbi_db=config["ncbi_db"]
    threads: 4
    conda:
        "../envs/malt.yaml"
    envmodules:
        *config["envmodules"]["malt"],
    message:
        "RUNNING MALT EXTRACT FOR SAMPLE {input.rma6}"
    shell:
        "time MaltExtract -i {input.rma6} -f def_anc -o {output.extract} --reads --threads {threads} --matches --minPI 85.0 --maxReadLength 0 --minComp 0.0 --meganSummary -r {params.ncbi_db} -t {input.node_list} -v"

#postprocessing.AMPS.r -m def_anc -r $OUT_DIR/${RMA6}_MaltExtract_output -t $THREADS -n $OUT_DIR/node_list.txt
rule Post_Processing:
    input:
        rma6="results/MALT/{sample}.trimmed.rma6",
        malt_extract_outdir=os.path.join("results/AUTHENTICATION/{sample}/{taxid}", "{sample}.trimmed.rma6_MaltExtract_output"),
        node_list="results/AUTHENTICATION/{sample}/{taxid}/node_list.txt",
    output:
        analysis="results/AUTHENTICATION/{sample}/{taxid}/{sample}.trimmed.rma6_MaltExtract_output/analysis.RData"
    threads: 4
    conda:
        "../envs/malt.yaml"
    envmodules:
        *config["envmodules"]["malt"],
    shell:
        "postprocessing.AMPS.r -m def_anc -r {input.malt_extract_outdir} -t {threads} -n {input.node_list}"
        

#head -2 $OUT_DIR/${RMA6}_MaltExtract_output/default/readDist/*.rma6_additionalNodeEntries.txt | tail -1 | cut -d ';' -f2 | sed 's/'_'/''/1' > $OUT_DIR/name.list
#REF_ID=$(cat $OUT_DIR/name.list)
#zgrep $REF_ID $IN_DIR/$SAM | uniq > $OUT_DIR/${REF_ID}.sam
#samtools view -bS $OUT_DIR/${REF_ID}.sam > $OUT_DIR/${REF_ID}.bam
#samtools sort $OUT_DIR/${REF_ID}.bam > $OUT_DIR/${REF_ID}.sorted.bam
#samtools index $OUT_DIR/${REF_ID}.sorted.bam
#samtools depth -a $OUT_DIR/${REF_ID}.sorted.bam > $OUT_DIR/${REF_ID}.breadth_of_coverage
#seqtk subseq $MALT_FASTA $OUT_DIR/name.list > $OUT_DIR/${REF_ID}.fasta
def get_ref_id(filename):
    f = open(filename)
    lines = f.readlines()
    ref_id = lines[1].split(";")[1][1:]
    
    return ref_id

rule Breadth_Of_Coverage:
    output:
        name_list="results/AUTHENTICATION/{sample}/{taxid,[0-9]+}/name.list",
        bam="results/AUTHENTICATION/{sample}/{taxid,[0-9]+}/{taxid}.sorted.bam",
    input:
        extract=os.path.join("results/AUTHENTICATION/{sample}/{taxid}", "{sample}.trimmed.rma6_MaltExtract_output"),
        dir="results/AUTHENTICATION/{sample}/{taxid}/",
        sam="results/MALT/{sample}.trimmed.sam.gz",
    params:
        malt_fasta=config["malt_nt_fasta"],
        #ref_id="ref_64"
    message:
        "COMPUTING BREADTH OF COVERAGE, EXTRACTING REFERENCE SEQUENCE FOR VISUALIZING ALIGNMENTS WITH IGV"
    conda:
        "../envs/malt.yaml"
    envmodules:
        *config["envmodules"]["malt"],
    shell:
        "head -2  results/AUTHENTICATION/{wildcards.sample}/{wildcards.taxid}/{wildcards.sample}.trimmed.rma6_MaltExtract_output/default/readDist/*.rma6_additionalNodeEntries.txt | tail -1 | cut -d ';' -f2 | sed 's/'_'/''/1' > {output.name_list};"
        "REF_ID=$(cat {output.name_list}); "
        "zgrep $REF_ID {input.sam} | uniq > {input.dir}/{wildcards.taxid}.sam; "
        "samtools view -bS {input.dir}/{wildcards.taxid}.sam > {input.dir}/{wildcards.taxid}.bam; "
        "samtools sort {input.dir}/{wildcards.taxid}.bam > {output.bam}; "
        "samtools index {input.dir}/{wildcards.taxid}.sorted.bam; "
        "samtools depth -a {input.dir}/{wildcards.taxid}.sorted.bam > {input.dir}/{wildcards.taxid}.breadth_of_coverage; "
        "seqtk subseq {params.malt_fasta} {output.name_list} > {input.dir}/{wildcards.taxid}.fasta"
        
#samtools view $OUT_DIR/${REF_ID}.sorted.bam | awk '{print length($10)}' > $OUT_DIR/${REF_ID}.read_length.txt
rule Read_Length_Distribution:
    output:
        distribution="results/AUTHENTICATION/{sample}/{taxid,[0-9]+}/{taxid}.read_length.txt"
    input:
        bam="results/AUTHENTICATION/{sample}/{taxid}/{taxid}.sorted.bam"
    message:
        "COMPUTING READ LENGTH DISTRIBUTION"
    conda:
        "../envs/malt.yaml"
    envmodules:
        *config["envmodules"]["malt"],
    shell:
        "samtools view {input.bam} | awk '{{ print length($10) }}' > {output.distribution}"

# samtools view -h $OUT_DIR/${REF_ID}.sorted.bam | pmdtools --printDS > $OUT_DIR/${REF_ID}.PMDscores.txt
rule PMD_scores:
    output:
        scores="results/AUTHENTICATION/{sample}/{taxid,[0-9]+}/{taxid}.PMDscores.txt"
    input:
        bam="results/AUTHENTICATION/{sample}/{taxid}/{taxid}.sorted.bam"
    message:
        "COMPUTING PMD SCORES"
    conda:
        "../envs/malt.yaml"
    envmodules:
        *config["envmodules"]["malt"],
    shell:
        "samtools view -h {input.bam} | pmdtools --printDS > {output.scores}"

# Rscript $AUTH_R_DIR/authentic.R $TAXID $IN_DIR $RMA6 $OUT_DIR
rule Authentication_Plots:
    output: 
        #dir=directory("results/AUTHENTICATION/{sample}/{taxid,[0-9]+}"),
        plot=os.path.join("results/AUTHENTICATION/{sample}/{taxid,[0-9]+}", "authentic_Sample_{sample}.trimmed.rma6_TaxID_{taxid}.pdf")
    input:
        rma6="results/MALT/{sample}.trimmed.rma6",
        node_list="results/AUTHENTICATION/{sample}/{taxid}/node_list.txt",
        name_list="results/AUTHENTICATION/{sample}/{taxid}/name.list",
        distribution="results/AUTHENTICATION/{sample}/{taxid}/{taxid}.read_length.txt",
        scores="results/AUTHENTICATION/{sample}/{taxid}/{taxid}.PMDscores.txt"
    params:
        exe=WORKFLOW_DIR / "scripts/authentic.R",
    message:
        "MAKING AUTHENTICATION AND VALIDATION PLOTS"
    conda:
        "../envs/malt.yaml"
    envmodules:
        *config["envmodules"]["malt"],
    shell:
        "Rscript {params.exe} {wildcards.taxid} results/MALT {wildcards.sample}.trimmed.rma6 results/AUTHENTICATION/{wildcards.sample}/{wildcards.taxid}"

rule Deamination:
    output:
        pmd=os.path.join("results/AUTHENTICATION/{sample}/{taxid,[0-9]+}", "PMD_temp.txt"),
    input:
        bam="results/AUTHENTICATION/{sample}/{taxid}/{taxid}.sorted.bam"
    message:
        "INFERRING DEAMINATION PATTERN FROM CPG SITES"
    conda:
        "../envs/malt.yaml"
    envmodules:
        *config["envmodules"]["malt"],
    shell:
        "samtools view {input.bam} | pmdtools --platypus > {output.pmd}; "
        "cd results/AUTHENTICATION/{wildcards.sample}/{wildcards.taxid}; "
        "R CMD BATCH $(which plotPMD); "
