#AUTHENTICATION WORKFLOW
#Example command line:
#scripts/./authentic.sh 632 input_dir input.rma6 input.sam.gz output_dir authentic.R_dir taxDB_dir ncbi_db malt_fasta

# See https://unix.stackexchange.com/questions/562018/run-every-bash-command-through-a-function-or-wrapper
function check {
    set -f # No double expand
    eval "$@" | tee >( tr -s "\n" " " )
    if [[ "${PIPESTATUS[0]}" != "0" ]]; then
        echo ERROR: $@
        exit 1
    fi
    return ${PIPESTATUS[0]}
}

arglen=$#
if [[ "$arglen" != '9' ]];
then
    echo "Wrong number of input arguments ($arglen)"
    echo "Usage: authentic.sh taxid indir rma6file samfile outdir scripts_dir taxdb ncbi_db malt_fasta"
    exit 0;
fi

TAXID=$1
IN_DIR=$2
RMA6=$(basename $3)
SAM=$(basename $4)
OUT_DIR=$5
AUTH_R_DIR=$6
TAXDB_DIR=$7
NCBI_DB=$8
MALT_FASTA=$9

#RUN MALT EXTRACT STATISTICS
echo "RUNNING MALT EXTRACT STATISTICS"
mkdir -p $OUT_DIR
check awk -v var="$TAXID" '{if($1==var)print$0}' $TAXDB_DIR/taxDB | cut -f3 > $OUT_DIR/node_list.txt
time check MaltExtract -i $IN_DIR/$RMA6 -f def_anc -o $OUT_DIR/${RMA6}_MaltExtract_output --reads --threads 4 --matches --minPI 85.0 --maxReadLength 0 --minComp 0.0 --meganSummary -r $NCBI_DB -t $OUT_DIR/node_list.txt -v
check postprocessing.AMPS.r -m def_anc -r $OUT_DIR/${RMA6}_MaltExtract_output -t 4 -n $OUT_DIR/node_list.txt

#COMPUTE BREADTH OF COVERAGE
echo "COMPUTING BREADTH OF COVERAGE"
check head -2 $OUT_DIR/${RMA6}_MaltExtract_output/default/readDist/*.rma6_additionalNodeEntries.txt | tail -1 | cut -d ';' -f2 | sed 's/'_'/''/1' > $OUT_DIR/name.list
REF_ID=$(cat $OUT_DIR/name.list)
check zgrep $REF_ID $IN_DIR/$SAM > $OUT_DIR/${REF_ID}.sam

check samtools view -bS $OUT_DIR/${REF_ID}.sam > $OUT_DIR/${REF_ID}.bam

check samtools sort $OUT_DIR/${REF_ID}.bam > $OUT_DIR/${REF_ID}.sorted.bam

check samtools index $OUT_DIR/${REF_ID}.sorted.bam

check samtools depth -a $OUT_DIR/${REF_ID}.sorted.bam > $OUT_DIR/${REF_ID}.breadth_of_coverage

#EXTRACT REFERENCE SEQUENCE FOR VISUALIZING ALIGNMENTS WITH IGV
echo "EXTRACTING REFERENCE SEQUENCE FOR VISUALIZING ALIGNMENTS WITH IGV"
check seqtk subseq $MALT_FASTA $OUT_DIR/name.list > $OUT_DIR/${REF_ID}.fasta

#COMPUTE READ LENGTH DISTRIBUTION
echo "COMPUTING READ LENGTH DISTRIBUTION"
check samtools view $OUT_DIR/${REF_ID}.sorted.bam | awk '{print length($10)}' > $OUT_DIR/${REF_ID}.read_length.txt

#COMPUTE PMD SCORES
echo "COMPUTING PMD SCORES"
check samtools view -h $OUT_DIR/${REF_ID}.sorted.bam | pmdtools --printDS > $OUT_DIR/${REF_ID}.PMDscores.txt

#MAKE AUTHENTICATION AND VALIDATION PLOTS
echo "MAKING AUTHENTICATION AND VALIDATION PLOTS"
check Rscript $AUTH_R_DIR/authentic.R $TAXID $IN_DIR $RMA6 $OUT_DIR

#INFER DEAMINATION PATTERN FROM CPG SITES (THAT ESCAPE USER TREATMENT)
echo "INFERRING DEAMINATION PATTERN FROM CPG SITES"
check samtools view $OUT_DIR/${REF_ID}.sorted.bam | pmdtools --platypus > $OUT_DIR/PMD_temp.txt

cd $OUT_DIR/
check R CMD BATCH $(which plotPMD)
