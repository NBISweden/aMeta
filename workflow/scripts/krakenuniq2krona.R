#This is a script for visualizing filtered KrakenUniq output with Krona. Run this script as:
#Rscipt krakenuniq2krona.R bal004.krakenuniq.output bal004.sequences.krakenuniq

args = commandArgs(trailingOnly=TRUE)

#Read and filter KrakenUniq output
df<-read.delim(args[1],comment.char="#",header=TRUE)
colnames(df)[1]<-"Pers_Reads"
print(paste0("Original data set dimensions: ",dim(df)[1]," and ",dim(df)[2]))
df<-df[df$kmers>1,]
print(paste0("Data set dimensions after breadth of coverage filter: ",dim(df)[1]," and ",dim(df)[2]))
df<-df[order(-df$Pers_Reads),]
df<-df[df$taxReads>1,]
df<-df[as.character(df$rank)=="species",]
print(paste0("Data set dimensions after depth of coverage filter: ",dim(df)[1]," and ",dim(df)[2]))
print(head(df))
write.table(df,file=paste0(args[1],".filtered"),col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")
write.table(as.character(df$taxID),file=paste0(args[1],"_taxIDs_kmers1000.txt"),col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")


#Select KrakenUniq classified reads corresponding to species remained after KrakenUniq output filtering
df<-read.delim(args[2],header=FALSE,sep="\t")
print(paste0("Original sequence data set dimensions: ",dim(df)[1]," and ",dim(df)[2]))
taxIDs<-scan(paste0(args[1],"_taxIDs_kmers1000.txt"),what="character")
df<-df[as.character(df$V3)%in%taxIDs,]
print(paste0("Sequence data set dimensions after selecting reads corresponding to filtered KrakenUniq output: ",dim(df)[1]," and ",dim(df)[2]))
write.table(df,file=paste0(args[2],"_kmers1000.txt"),col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")
