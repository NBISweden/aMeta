#This is a script for combining KrakenUniq outputs from multiple samples and producing species abundance matrix.
#Run this script as:
#Rscipt krakenuniq_abundance_matrix.R input_dir output_dir

args = commandArgs(trailingOnly=TRUE)

#Reading and merging krakenuniq outputs from multiple samples
krakenuniq_outputs<-list.files(path=args[1])
df<-list()
for(i in 1:length(krakenuniq_outputs))
{
 df[[i]]<-read.delim(paste0(args[1],"/",krakenuniq_outputs[i],"/krakenuniq.output"),comment.char="#",header=TRUE)
 names(df[[i]])[1]<-"Pers_Reads"
 df[[i]]$SAMPLE<-krakenuniq_outputs[i]
 df[[i]]<-na.omit(df[[i]])
 df[[i]]<-df[[i]][df[[i]]$kmers>1,]
 df[[i]]<-df[[i]][as.character(df[[i]]$rank)=="species",]
 df[[i]]<-df[[i]][df[[i]]$taxReads>1,]
}
merged<-Reduce(rbind,df)
merged$taxName<-trimws(as.character(merged$taxName))
print(head(merged,20))

unique_samples<-unique(merged$SAMPLE)
unique_species<-unique(merged$taxName)
unique_taxids<-unique(merged$taxID)

#Computing abundance matrix
abundance_matrix<-matrix(NA,ncol=length(unique_samples),nrow=length(unique_species))
for(i in 1:length(unique_species))
{
 for(j in 1:length(unique_samples))
 {
  if(length(merged[merged$taxName==unique_species[i] & merged$SAMPLE==unique_samples[j],]$taxReads)==0)
  {
   abundance_matrix[i,j]<-0
  }
  else
  {
   abundance_matrix[i,j]<-merged[merged$taxName==unique_species[i] & merged$SAMPLE==unique_samples[j],]$taxReads
  }
 }
}
rownames(abundance_matrix)<-unique_species
colnames(abundance_matrix)<-unique_samples
print(head(abundance_matrix))

system(paste0("mkdir ",args[2]))
write.table(abundance_matrix,file=paste0(args[2],"/krakenuniq_abundance_matrix.txt"),col.names=TRUE,row.names=TRUE,quote=FALSE,sep="\t")
write.table(unique_taxids,file=paste0(args[2],"/unique_species_taxid_list.txt"),col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")
write.table(unique_species,file=paste0(args[2],"/unique_species_names_list.txt"),col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")
