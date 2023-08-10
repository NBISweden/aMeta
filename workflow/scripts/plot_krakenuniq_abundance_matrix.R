#This is a script for plotting KrakenUniq abundance matrix.
#Run this script as:
#Rscript plot_krakenuniq_abundance_matrix.R in_dir out_dir

args = commandArgs(trailingOnly=TRUE)
in_dir<-as.character(args[1])
out_dir<-as.character(args[2])

library("pheatmap")
ku_abundance<-read.delim(paste0(in_dir,"/krakenuniq_abundance_matrix.txt"),header=TRUE,row.names=1,check.names=FALSE,sep="\t")


#FUNCTION FOR TUNING FONTSIZE ON MICROBIAL ABUNDANCE HEATMAP
my_fontsize<-function(ku_abundance)
{
  if(dim(ku_abundance)[1]<50){return(12)}
  else if(dim(ku_abundance)[1]>=50 & dim(ku_abundance)[1]<100){return(10)}
  else if(dim(ku_abundance)[1]>=100 & dim(ku_abundance)[1]<150){return(8)}
  else{return(6)}
}


#ABSOLUTE ABUNDANCE HEATMAP
pdf(paste0(out_dir,"/krakenuniq_absolute_abundance_heatmap.pdf"),paper="a4r",width=297,height=210)
if(dim(ku_abundance)[1]>1 & dim(ku_abundance)[2]>1)
{
  pheatmap(ku_abundance, display_numbers=TRUE,fontsize=my_fontsize(ku_abundance),
           main="KrakenUniq Absolute Microbial Abundance",cluster_rows=FALSE,cluster_cols=FALSE,number_format="%i")
}else
{
  pheatmap(ku_abundance, display_numbers=TRUE,fontsize=8,
           main="KrakenUniq Absolute Microbial Abundance",cluster_rows=FALSE,cluster_cols=FALSE,number_format="%i",breaks=c(0,1))
}
dev.off()


#NORMALIZE BY SEQUENCING SEPTH
for(i in 1:dim(ku_abundance)[2])
{
  ku_abundance[,i]<-ku_abundance[,i]/sum(ku_abundance[,i])
}


#ABSOLUTE ABUNDANCE HEATMAP
pdf(paste0(out_dir,"/krakenuniq_normalized_abundance_heatmap.pdf"),paper="a4r",width=297,height=210)
if(dim(ku_abundance)[1]>1 & dim(ku_abundance)[2]>1)
{
  pheatmap(ku_abundance, display_numbers=TRUE,fontsize=my_fontsize(ku_abundance),
           main="KrakenUniq Normalized Microbial Abundance",cluster_rows=FALSE,cluster_cols=FALSE,number_format="%.3f")
}else
{
  pheatmap(ku_abundance, display_numbers=TRUE,fontsize=8,
           main="KrakenUniq Normalized Microbial Abundance",cluster_rows=FALSE,cluster_cols=FALSE,number_format="%.3f",breaks=c(0,1))
}
dev.off()
