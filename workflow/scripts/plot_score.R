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
  taxa<-list.files(file.path(in_dir, i))
  if(length(taxa)!=0)
  {
    score<-list()
    for(j in taxa)
    {
      if(file.exists(file.path(in_dir, i, j, "authentication_scores.txt"))==TRUE)
      {
        score[[j]]<-read.delim(file.path(in_dir, i, j, "authentication_scores.txt"),header=TRUE,sep="\t")
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
score_matrix$ORGANISM<-NULL
score_matrix[is.na(score_matrix)]<-0
score_matrix
write.table(score_matrix,file=file.path(out_dir, "overview_heatmap_scores.txt"),col.names=TRUE,row.names=TRUE,quote=FALSE,sep="\t")
                     

library("pheatmap")
pdf(file.path(out_dir, "overview_heatmap_scores.pdf"),paper="a4r",width=297,height=210)
if(dim(score_matrix)[1]>1 & dim(score_matrix)[2]>1)
{
  pheatmap(score_matrix, display_numbers=TRUE,fontsize=8,main="Ancient microbiome profiling overview",
           cluster_rows=FALSE,cluster_cols=FALSE,number_format="%i")
}else
{
  pheatmap(score_matrix, display_numbers=TRUE,fontsize=8,main="Ancient microbiome profiling overview",
           cluster_rows=FALSE,cluster_cols=FALSE,number_format="%i",breaks=c(0,1))
}
dev.off()
