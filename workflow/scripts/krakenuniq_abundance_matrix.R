#This is a script for combining KrakenUniq outputs from multiple samples and producing species abundance matrix.
#Run this script as:
#Rscipt krakenuniq_abundance_matrix.R input_dir output_dir n_unique_kmers n_tax_reads

args = commandArgs(trailingOnly=TRUE)
input_dir<-as.character(args[1])
output_dir<-as.character(args[2])
n_unique_kmers<-as.integer(args[3])
n_tax_reads<-as.integer(args[4])

#Reading and merging krakenuniq outputs from multiple samples
krakenuniq_outputs<-list.files(path=input_dir)
df<-list()
for(i in 1:length(krakenuniq_outputs))
{
 df[[i]]<-read.delim(paste0(input_dir,"/",krakenuniq_outputs[i],"/krakenuniq.output.filtered"),comment.char="#",header=TRUE)
 if(dim(df[[i]])[1]!=0)
 {
  df[[i]]$SAMPLE<-krakenuniq_outputs[i]
  df[[i]]<-na.omit(df[[i]])
 }
 else{next}
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
abundance_matrix<-abundance_matrix[order(rownames(abundance_matrix)),]
print(head(abundance_matrix))

system(paste0("mkdir ",output_dir))
write.table(abundance_matrix,file=paste0(output_dir,"/krakenuniq_abundance_matrix.txt"),col.names=TRUE,row.names=TRUE,quote=FALSE,sep="\t")
write.table(unique_taxids,file=paste0(output_dir,"/unique_species_taxid_list.txt"),col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")
write.table(unique_species,file=paste0(output_dir,"/unique_species_names_list.txt"),col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")
