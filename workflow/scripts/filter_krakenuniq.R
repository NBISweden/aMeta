#This is a script for filtering KrakenUniq output. Run this script as:
#Rscipt filter_krakenuniq.R krakenuniq.output

args = commandArgs(trailingOnly=TRUE)

#Reading KrakenUniq output
df <- read.delim(args[1],comment.char="#",header=TRUE)

#Filtering KrakenUniq output with respect to depth (taxReads) and breadth (kmers) of coverage
colnames(df)[1] <- "Pers_Reads"
df <- df[df$kmers>1,]
df <- df[as.character(df$rank) == "species",]
df <- df[df$taxReads>1,]
write.table(df, file = paste0(args[1],".filtered"), col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

#Overlapping the filtered KrakenUniq output with the list of pathogens
pathogens <- read.delim(args[2], header=FALSE)
df <- df[as.character(df$taxID)%in%pathogens$V1,]
write.table(df, file = paste0(args[1],".pathogens"), col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
