#This is a script for filtering KrakenUniq output. Run this script as:
#Rscipt filter_krakenuniq.R krakenuniq_output pathogens_list n_unique_kmers n_tax_reads

args = commandArgs(trailingOnly=TRUE)
krakenuniq_output<-as.character(args[1])
pathogens_list<-as.character(args[2])
n_unique_kmers<-as.integer(args[3])
n_tax_reads<-as.integer(args[4])


#Reading KrakenUniq output
df <- read.delim(krakenuniq_output,comment.char="#",header=TRUE)

#Filtering KrakenUniq output with respect to depth (taxReads) and breadth (kmers) of coverage
colnames(df)[1] <- "Pers_Reads"
df <- df[df$kmers>n_unique_kmers,]
df <- df[as.character(df$rank) == "species",]
df <- df[df$taxReads>n_tax_reads,]
write.table(df, file = paste0(krakenuniq_output,".filtered"), col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

#Overlapping the filtered KrakenUniq output with the list of pathogens
pathogens <- read.delim(pathogens_list, header=FALSE)
df <- df[as.character(df$taxID)%in%pathogens$V1,]
write.table(df, file = paste0(krakenuniq_output,".pathogens"), col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
