#This is a script for plotting a heatmap overview of aMeta authentication scores.
#Run this script as:
#Rscript plot_score.R in_dir out_dir

#in_dir<-"aMeta/results/AUTHENTICATION"
#out_dir<-"aMeta/results"

args = commandArgs(trailingOnly=TRUE)
in_dir<-as.character(args[1])
out_dir<-as.character(args[2])

samples<-list.files(in_dir)

score_per_sample<-list()
for(i in samples)
{
  taxa<-list.files(paste0(in_dir,"/",i))
  if(length(taxa)!=0)
  {
    score<-list()
    for(j in taxa)
    {
      if(file.exists(paste0(in_dir,"/",i,"/",j,"/authentication_scores.txt"))==TRUE)
      {
        score[[j]]<-read.delim(paste0(in_dir,"/",i,"/",j,"/authentication_scores.txt"),header=TRUE,sep="\t")
      }else
      {
        next
      }
      
    }
    score_per_sample[[i]]<-Reduce(rbind,score)
    names(score_per_sample[[i]])[2]<-i
  }else
  {
    next
  }
}
score_matrix<-Reduce(function(x,y) merge(x,y,all=TRUE,by="ORGANISM"),score_per_sample)
rownames(score_matrix)<-score_matrix$ORGANISM
score_matrix<-score_matrix[order(rownames(score_matrix)),]
score_matrix$ORGANISM<-NULL
score_matrix[is.na(score_matrix)]<-0
score_matrix<-score_matrix[grepl("Homo sapiens",rownames(score_matrix))==FALSE,]
score_matrix
write.table(score_matrix,file=paste0(out_dir,"/overview_heatmap_scores.txt"),col.names=TRUE,row.names=TRUE,quote=FALSE,sep="\t")


#FUNCTION FOR TUNING FONTSIZE ON SCORE HEATMAP
my_fontsize<-function(score_matrix)
{
  if(dim(score_matrix)[1]<50){return(12)}
  else if(dim(score_matrix)[1]>=50 & dim(score_matrix)[1]<100){return(10)}
  else if(dim(score_matrix)[1]>=100 & dim(score_matrix)[1]<150){return(8)}
  else{return(6)}
}

library("pheatmap")
pdf(paste0(out_dir,"/overview_heatmap_scores.pdf"),paper="a4r",width=297,height=210)
if(dim(score_matrix)[1]>1 & dim(score_matrix)[2]>1)
{
  pheatmap(score_matrix, display_numbers=TRUE,fontsize=my_fontsize(score_matrix),main="Ancient microbiome profiling overview",
           cluster_rows=FALSE,cluster_cols=FALSE,number_format="%i")
}else
{
  pheatmap(score_matrix, display_numbers=TRUE,fontsize=12,main="Ancient microbiome profiling overview",
           cluster_rows=FALSE,cluster_cols=FALSE,number_format="%i",breaks=c(0,1))
}
dev.off()

