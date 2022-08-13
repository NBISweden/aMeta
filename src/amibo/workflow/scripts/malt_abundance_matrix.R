#This is a script for combining Malt microbial quantifications from multiple samples and producing species abundance matrix.
#Run this script as:
#Rscipt malt_abundance_matrix.R input_dir output_dir

args = commandArgs(trailingOnly=TRUE)

#Read species scientific names to be added as rownames to the abundance matrix
species_names<-read.delim("results/KRAKENUNIQ_ABUNDANCE_MATRIX/unique_species_names_list.txt",header=FALSE,sep="\t")

#Read Malt microbial quantifications (each file is one vector of counts) and merge them horizontally
files<-list.files(path=args[1])
df<-list()
for(i in 1:length(files))
{
df[[i]]<-as.numeric(scan(paste0(args[1],"/",files[i],"/sam_counts.txt"),what="character"))
}
merged<-Reduce(cbind,df)
colnames(merged)<-gsub(".sam_counts","",files)
rownames(merged)<-as.character(species_names$V1)

#Write the resulting Malt abundance matrix to file
write.table(merged,file=paste0(args[2],"/malt_abundance_matrix_sam.txt"),col.names=TRUE,row.names=TRUE,quote=FALSE,sep="\t")
