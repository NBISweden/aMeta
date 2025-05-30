#MAKING VALIDATION PLOTS

#ARGUMENTS OF SCRIPT
args = commandArgs(trailingOnly=TRUE)
taxid<-args[1]
RMA6<-args[2]
out_dir<-args[3]

options(warn = - 1)

pdf(paste0(out_dir,"/","authentic_Sample_",RMA6,"_TaxID_",taxid,".pdf"),paper="a4r",width=297,height=210)
par(mfrow=c(3,3))

#HELP INFORMATION
organism<-readLines(paste0(out_dir,"/node_list.txt")) #scientific name of oranism, extracted automatically from NCBI NT by taxID
MaltExtract_output_path<-paste0(out_dir,"/MaltExtract_output") #path to MaltExtract output directory

#EDIT DISTANCE FOR ALL READS
df<-read.delim(paste0(MaltExtract_output_path,"/default/editDistance/",RMA6,"_editDistance.txt"),header=TRUE,check.names=FALSE,row.names=1,sep="\t")
df$higher<-NULL
barplot(as.numeric(df[1,]),names=c(0:10),xlab="Number of mismatches",ylab="Number of reads",main="Edit distance: all reads")

#EDIT DISTANCE FOR ANCIENT READS
df<-read.delim(paste0(MaltExtract_output_path,"/ancient/editDistance/",RMA6,"_editDistance.txt"),header=TRUE,check.names=FALSE,row.names=1,sep="\t")
df$higher<-NULL
barplot(as.numeric(df[1,]),names=c(0:10),xlab="Number of mismatches",ylab="Number of reads",main="Edit distance: ancient reads")

#BREADTH OF COVERAGE
# FIXME: empty due to failed samtools sort
df<-read.delim(paste0(out_dir,"/breadth_of_coverage"),header=FALSE,sep="\t")
N_tiles<-100
step=(max(df$V2)-min(df$V2))/N_tiles
tiles<-c(0:N_tiles)*step
V4<-vector()
for(i in 1:length(tiles))
{
  df_temp<-df[df$V2>=tiles[i] & df$V2<tiles[i+1],]
  V4<-append(V4,rep(sum(df_temp$V3>0)/length(df_temp$V3),dim(df_temp)[1]))
}
V4[is.na(V4)]<-0
df$V4<-V4
plot(df$V4~df$V2,type="s",xlab="Genome position",ylab="Fraction of covered genome",main=paste0("Evenness of coverage: ",readLines(paste0(out_dir,"/name_list.txt"))," reference"))
abline(h=0,col="red",lty=2)
mtext(paste0("Breadth of coverage: ",round((sum(df$V3>0)/length(df$V3))*100,2),"% of genome covered"),cex=0.8)

#DAMAGE PATTERN
dam <- read.table(paste0(MaltExtract_output_path,"/default/damageMismatch/",RMA6,"_damageMismatch.txt"),header=T,row.names=1,check.names=F,stringsAsFactors=F,comment.char='')
dam <- dam[1, ]
maxY <- (max(dam[, -dim(dam)[2] ]) * 1.1)
plot("",col="grey",ylim=c(0,maxY), xlim=c(1,20),type="l",lwd=2,ylab="C2T / G2A transition rate",xlab="Read Position",xaxt='n',main="Deamination pattern plot")
lines(x=1:10,dam[1,21:30],col='grey',lwd=1.5)
lines(x=11:20,dam[1,31:40],col='grey',lwd=1.5)
lines(x=1:10,y=dam[1,1:10],col='red',lwd=1.5)
lines(x=11:20,dam[1,11:20],col='blue',lwd=1.5)
axis(1, at=seq(1,20,2), labels=c(seq(1,10,2),seq(-10,-1,2)))

#READ LENGTH DISTRIBUTION
df<-as.numeric(scan(paste0(out_dir,"/read_length.txt"),what="character"))
hist(df,breaks=length(df),xlab="Read length",main="Read length distribution")

#PMD SCORES DISTRIBUTION
if(file.info(paste0(out_dir,"/PMDscores.txt"))$size!=0)
{
df<-read.delim(paste0(out_dir,"/PMDscores.txt"),header=FALSE,sep="\t")
hist(df$V4,breaks=length(df$V4),main="Histogram of PMD scores",xlab="PMDscores")
abline(v=3,col="red",lty=2)
}else{
plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', main="Histogram of PMD scores")
text(x = 0.5, y = 0.5,paste("Not enough reads"),cex=1.6)
}

#PERCENT IDENTITY
df<-read.delim(paste0(MaltExtract_output_path,"/default/percentIdentity/",RMA6,"_percentIdentity.txt"),header=TRUE,check.names=FALSE,row.names=1,sep="\t")
barplot(as.numeric(df[1,]),names=colnames(df),xlab="Percent identity",ylab="Number of reads",main="Reads mapped with identity to reference")
mtext(paste0("Average nucleotide identity (ANI) = ",round((df[1,"90"]*90+df[1,"95"]*95+df[1,"100"]*100)/(df[1,"90"]+df[1,"95"]+df[1,"100"]),1),"%"),cex=0.8)

#TABLE OF TOP MAPPING REFERENCES
library("gridBase")
library("gridExtra")
ar<-read.table(paste0(MaltExtract_output_path,"/default/readDist/",RMA6,"_additionalNodeEntries.txt"),header=T,row.names=1,check.names=F,stringsAsFactors=F,comment.char='',sep="\t",quote="")
ar<-ar[1, ]
ar<-paste(sub(";_TOPREFPERCREADS"," ",ar) ,"%",sep="")
plot.new()
mytheme<-gridExtra::ttheme_default(base_size=8) #https://github.com/baptiste/gridextra/wiki/tableGrob
grid.table(ar, vp=gridBase::baseViewports()$figure,theme=mytheme)

#TABLE OF USEFUL STATISTICS
library("gridBase")
library("gridExtra")
rd<-read.table(paste0(MaltExtract_output_path,"/default/readDist/",RMA6,"_alignmentDist.txt"),header=T,row.names=1,check.names=F,stringsAsFactors=F,comment.char='')
rd<-rd[1, ]
topNode<-rownames(rd)
dam<-read.table(paste0(MaltExtract_output_path,"/default/damageMismatch/",RMA6,"_damageMismatch.txt"),header=T,row.names=1,check.names=F,stringsAsFactors=F,comment.char='')
idxN<-dim(dam)[2]
dam<-dam[1, ]
ld<-read.table(paste0(MaltExtract_output_path,"/default/readDist/",RMA6,"_readLengthStat.txt"),header=T,row.names=1,check.names=F,stringsAsFactors=F,comment.char='')
ld<-paste( round(ld[1, 'Mean' ],0),' (',round(ld[1, 'StandardDev' ],3),')',sep="")
ds<-read.table(paste0(MaltExtract_output_path,"/default/filterInformation/",RMA6,"_filterTable.txt"),header=T,row.names=1,check.names=F,stringsAsFactors=F,comment.char='')
ds<-ds[1, "turnedOn?"]
data<-c(topNode,ar[1],rd[topNode,'TotalAlignmentsOnReference'],rd[topNode,'nonDuplicatesonReference'],rd[topNode,'uniquePerReference'],rd[topNode,'nonStacked'],ds,round(dam[topNode,'C>T_1'],4),round(dam[topNode,'G>A_20'],4),ld)
data<-cbind( c('Node','Top Reference','all reads','nonDup','readDis','nonStacked','destacking?','C>T_1','G>A_-1','mean length (sd)'), data)
colnames(data)=NULL; rownames(data)=NULL
plot.new()
mytheme<-gridExtra::ttheme_default(base_size=8) #https://github.com/baptiste/gridextra/wiki/tableGrob
grid.table(data, vp=gridBase::baseViewports()$figure,theme=mytheme)

mtext(paste0("Organism: ",organism,", taxID: ",taxid,", rma6-file: ",RMA6),outer=TRUE,cex=0.8,line=-1)
dev.off()


system(paste0("convert -density 300 -background white -alpha remove ",paste0(out_dir,"/","authentic_Sample_",RMA6,"_TaxID_",taxid,".pdf")," ",paste0(out_dir,"/","authentic_Sample_",RMA6,"_TaxID_",taxid,".png")))
