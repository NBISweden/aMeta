#SCRIPT FOR COMPUTING ANCIENT_METAGENOME SCORE PER MICROBE
#RUN SCRIPT AS: Rscript score.R RMA6_FILE_NAME OUT_DIR DIR_NAME_LIST

args<-commandArgs(trailingOnly=TRUE)
RMA6<-args[1]
out_dir<-args[2]
dir_name_list<-args[3]
options(warn = - 1)

#out_dir<-"/home/nikolay/WABI/A_Gotherstrom/Manuscript/Method_Paper/HOPS_vs_AncientMetagenome/632"
#RMA6<-"simulation_s1.trimmed.rma6"

organism<-suppressWarnings(readLines(paste0(out_dir,"/node_list.txt"))) #scientific name of oranism, extracted automatically from NCBI NT by taxID
RefID<-suppressWarnings(readLines(paste0(dir_name_list,"/name.list"))) #sequence ID of reference sequence that has most of reads mapped to it
MaltExtract_output_path<-paste0(out_dir,"/",RMA6,"_MaltExtract_output") #path to MaltExtract output directory
rd<-suppressWarnings(read.table(paste0(MaltExtract_output_path,"/default/readDist/",RMA6,"_alignmentDist.txt"),header=T,row.names=1,check.names=F,stringsAsFactors=F,comment.char=''))
rd<-rd[1, ]
topNode<-rownames(rd)

total_score<-0

#DAMAGE PATTERN
dam<-suppressWarnings(read.table(paste0(MaltExtract_output_path,"/default/damageMismatch/",RMA6,"_damageMismatch.txt"),header=T,row.names=1,check.names=F,stringsAsFactors=F,comment.char=''))
dam<-dam[1, ]
if(round(dam[topNode,'C>T_1'],4)>0.05){total_score<-total_score+1}
if(round(dam[topNode,'G>A_20'],4)>0.05){total_score<-total_score+1}


#EVENNES OF COVERAGE
df<-try(read.delim(paste0(out_dir,"/",RefID,".breadth_of_coverage"),header=FALSE,sep="\t"),silent=TRUE)
if(!inherits(df,'try-error')){
N_tiles<-100
step=(max(df$V2,na.rm=TRUE)-min(df$V2,na.rm=TRUE))/N_tiles
tiles<-c(0:N_tiles)*step
V4<-vector()
for(i in 1:length(tiles))
{
  df_temp<-df[df$V2>=tiles[i] & df$V2<tiles[i+1],]
  V4<-append(V4,rep(sum(df_temp$V3>0)/length(df_temp$V3),dim(df_temp)[1]))
}
V4[is.na(V4)]<-0
if(sum(V4<0.01,na.rm=TRUE)/step<3){total_score<-total_score+3}
}else{
total_score<-total_score+3
}


#EDIT DISTANCE FOR ALL READS
df<-suppressWarnings(read.delim(paste0(MaltExtract_output_path,"/default/editDistance/",RMA6,"_editDistance.txt"),header=TRUE,check.names=FALSE,row.names=1,sep="\t"))
df$higher<-NULL
if(as.numeric(df[1,1])>as.numeric(df[1,2]) & as.numeric(df[1,2])>as.numeric(df[1,3]) & as.numeric(df[1,3])>as.numeric(df[1,4]) & as.numeric(df[1,4])>as.numeric(df[1,5])){total_score<-total_score+1}

#EDIT DISTANCE FOR ANCIENT READS
df<-suppressWarnings(read.delim(paste0(MaltExtract_output_path,"/ancient/editDistance/",RMA6,"_editDistance.txt"),header=TRUE,check.names=FALSE,row.names=1,sep="\t"))
df$higher<-NULL
if(as.numeric(df[1,2])>as.numeric(df[1,3]) & as.numeric(df[1,3])>as.numeric(df[1,4]) & as.numeric(df[1,4])>as.numeric(df[1,5])){total_score<-total_score+1}


#PMD SCORES
df<-try(read.delim(paste0(out_dir,"/",RefID,".PMDscores.txt"),header=FALSE,sep="\t"),silent=TRUE)
if(!inherits(df,'try-error')){
if(sum(df$V4>3)/dim(df)[1]>0.1){total_score<-total_score+1}
}else{
total_score<-total_score+1
}


#READ LENGTH DISTRIBUTION
if(file.exists(paste0(out_dir,"/",RefID,".read_length.txt")))
{
df<-as.numeric(readLines(paste0(out_dir,"/",RefID,".read_length.txt")))
if(length(df)>0){
if(sum(df<100)/length(df)>0.9){total_score<-total_score+1}
}else{
total_score<-total_score+1
}
}else{
total_score<-total_score+1
}


#DEPTH OF COVERAGE
if(rd[topNode,'TotalAlignmentsOnReference']>200){total_score<-total_score+1}


print(paste0(organism," score: ",total_score))



