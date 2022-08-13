#This is a script for visualizing filtered KrakenUniq output with Krona. Run this script as:
#Rscikrakenuniq_outputpt krakenuniq2krona.R krakenuniq_output sequences_krakenuniq n_unique_kmers n_tax_reads

args = commandArgs(trailingOnly=TRUE)
krakenuniq_output<-as.character(args[1])
sequences_krakenuniq<-as.character(args[2])
n_unique_kmers<-as.integer(args[3])
n_tax_reads<-as.integer(args[4])

#Read and filter KrakenUniq output
df<-read.delim(krakenuniq_output,comment.char="#",header=TRUE)
colnames(df)[1]<-"Pers_Reads"
print(paste0("Original data set dimensions: ",dim(df)[1]," and ",dim(df)[2]))
df<-df[df$kmers>n_unique_kmers,]
print(paste0("Data set dimensions after breadth of coverage filter: ",dim(df)[1]," and ",dim(df)[2]))
df<-df[order(-df$Pers_Reads),]
df<-df[df$taxReads>n_tax_reads,]
df<-df[as.character(df$rank)=="species",]
print(paste0("Data set dimensions after depth of coverage filter: ",dim(df)[1]," and ",dim(df)[2]))
print(head(df))
write.table(df,file=paste0(krakenuniq_output,".filtered"),col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")
write.table(as.character(df$taxID),file=paste0(krakenuniq_output,"_taxIDs_kmers1000.txt"),col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")


#Select KrakenUniq classified reads corresponding to species remained after KrakenUniq output filtering
df<-read.delim(sequences_krakenuniq,header=FALSE,sep="\t")
print(paste0("Original sequence data set dimensions: ",dim(df)[1]," and ",dim(df)[2]))
taxIDs<-scan(paste0(krakenuniq_output,"_taxIDs_kmers1000.txt"),what="character")
df<-df[as.character(df$V3)%in%taxIDs,]
print(paste0("Sequence data set dimensions after selecting reads corresponding to filtered KrakenUniq output: ",dim(df)[1]," and ",dim(df)[2]))
write.table(df,file=paste0(sequences_krakenuniq,"_kmers1000.txt"),col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")
