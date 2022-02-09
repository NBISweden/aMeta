#AUTHENTICATION WORKFLOW
#Example command line:
#scripts/./authentic.sh 632 input_dir input.rma6 input.sam.gz output_dir authentic.R_dir
#scripts/./authentic.sh 632 test_in Gok2c_GCAAGAT_L004_merged_001.130620_SN344_0250_AD2538ACXX.rma6 Gok2c_GCAAGAT_L004_merged_001.130620_SN344_0250_AD2538ACXX.sam.gz test_out scripts

TAXID=$1
IN_DIR=$2
RMA6=$(basename $3)
SAM=$(basename $4)
OUT_DIR=$5
AUTH_R_DIR=$6

ml bioinfo-tools seqtk samtools java python/2.7.6 R/3.6.0 R_packages/3.6.0

#RUN MALT EXTRACT STATISTICS
echo "RUNNING MALT EXTRACT STATISTICS"
mkdir $OUT_DIR
awk -v var="$TAXID" '{if($1==var)print$0}' /crex/proj/uppstore2018095/private/NBIS_Demo/DBDIR_KrakenUniq_Full_NT/taxDB.orig | cut -f3 > $OUT_DIR/node_list.txt
time MaltExtract -i $IN_DIR/$RMA6 -f def_anc -o $OUT_DIR/${RMA6}_MaltExtract_output --reads --threads 4 --matches --minPI 85.0 --maxReadLength 0 --minComp 0.0 --meganSummary -r /crex/proj/uppstore2018095/private/NBIS_Demo/HOPS/HOPS_github_clone/Resources -t $OUT_DIR/node_list.txt -v
postprocessing.AMPS.r -m def_anc -r $OUT_DIR/${RMA6}_MaltExtract_output -t 4 -n $OUT_DIR/node_list.txt

#COMPUTE BREADTH OF COVERAGE
echo "COMPUTING BREADTH OF COVERAGE"
head -2 $OUT_DIR/${RMA6}_MaltExtract_output/default/readDist/*.rma6_additionalNodeEntries.txt | tail -1 | cut -d ';' -f2 | sed 's/_//g' > $OUT_DIR/name.list
REF_ID=$(cat $OUT_DIR/name.list)
zgrep $REF_ID $IN_DIR/$SAM > $OUT_DIR/${REF_ID}.sam
samtools view -bS $OUT_DIR/${REF_ID}.sam > $OUT_DIR/${REF_ID}.bam
samtools sort $OUT_DIR/${REF_ID}.bam > $OUT_DIR/${REF_ID}.sorted.bam
samtools index $OUT_DIR/${REF_ID}.sorted.bam
samtools depth -a $OUT_DIR/${REF_ID}.sorted.bam > $OUT_DIR/${REF_ID}.breadth_of_coverage

#EXTRACT REFERENCE SEQUENCE FOR VISUALIZING ALIGNMENTS WITH IGV
echo "EXTRACTING REFERENCE SEQUENCE FOR VISUALIZING ALIGNMENTS WITH IGV"
seqtk subseq /proj/snic2018-8-150/uppstore2018095/private/NBIS_Demo/MALT_Full_NT_Database/20210127-003704_nt/library/nt/library.fna $OUT_DIR/name.list > $OUT_DIR/${REF_ID}.fasta

#COMPUTE READ LENGTH DISTRIBUTION
echo "COMPUTING READ LENGTH DISTRIBUTION"
samtools view $OUT_DIR/${REF_ID}.sorted.bam | awk '{print length($10)}' > $OUT_DIR/${REF_ID}.read_length.txt

#COMPUTE PMD SCORES
echo "COMPUTING PMD SCORES"
samtools view -h $OUT_DIR/${REF_ID}.sorted.bam | python2 /crex/proj/uppstore2018095/private/NBIS_Demo/workshop/PMDtools/pmdtools.0.60.py --printDS > $OUT_DIR/${REF_ID}.PMDscores.txt

#MAKE AUTHENTICATION AND VALIDATION PLOTS
echo "MAKING AUTHENTICATION AND VALIDATION PLOTS"
Rscript $AUTH_R_DIR/authentic.R $TAXID $IN_DIR $RMA6 $OUT_DIR

#INFER DEAMINATION PATTERN FROM CPG SITES (THAT ESCAPE USER TREATMENT)
echo "INFERRING DEAMINATION PATTERN FROM CPG SITES"
samtools view $OUT_DIR/${REF_ID}.sorted.bam | python2 /crex/proj/uppstore2018095/private/NBIS_Demo/workshop/PMDtools/pmdtools.0.60.py --platypus --requirebaseq 30 > $OUT_DIR/PMD_temp.txt
cd $OUT_DIR/
R CMD BATCH /crex/proj/uppstore2018095/private/NBIS_Demo/workshop/PMDtools/plotPMD.v2.R
