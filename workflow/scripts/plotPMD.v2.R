###### Companion script for plotting the results of PMDtools (Skoglund, Northoff, Shunkov, Dervianko, Paabo, Krause, Jakobsson, 2014, PNAS)
###### Contact: pontus.skoglund@gmail.com 
###### 
###### Usage example: 
######        samtools view mybam.bam | python pmdtools.py --deamination > PMD_temp.txt
######        R CMD BATCH plotPMD.R
######        cp PMD_plot.pdf mybam_PMD_plot.pdf



pdf("PMD_plot.frag.pdf",width=12,height=4,useDingbats=FALSE)


par(mfrow=c(1,3))


data<-read.table("PMD_temp.txt",header= TRUE)



read_position <-data$z

plot(read_position,read_position,xlab="Distance from 5' end of sequence read",ylab="Mismatch frequency",cex.axis=1,cex.lab=1.5,cex=3,type="n",ylim=c(0,0.50),xlim=c(0,max(read_position)))
lines(read_position,data$CG5,col="black",lwd=2)
lines(read_position,data$CA5,col="black",lwd=2)
lines(read_position,data$GT5,col="black",lwd=2)
lines(read_position,data$GC5,col="black",lwd=2)
lines(read_position,data$AG5,col="black",lwd=2)
lines(read_position,data$AT5,col="black",lwd=2)
lines(read_position,data$AC5,col="black",lwd=2)
lines(read_position,data$TC5,col="black",lwd=2)
lines(read_position,data$TG5,col="black",lwd=2)
lines(read_position,data$TA5,col="black",lwd=2)

lines(read_position,data$CG_CpG_5,col="black",lwd=2,lty=2)
lines(read_position,data$CA_CpG_5,col="black",lwd=2,lty=2)
lines(read_position,data$GT_CpG_5,col="black",lwd=2,lty=2)
lines(read_position,data$GC_CpG_5,col="black",lwd=2,lty=2)
lines(read_position,data$AG_CpG_5,col="black",lwd=2,lty=2)
lines(read_position,data$AT_CpG_5,col="black",lwd=2,lty=2)
lines(read_position,data$AC_CpG_5,col="black",lwd=2,lty=2)
lines(read_position,data$TC_CpG_5,col="black",lwd=2,lty=2)
lines(read_position,data$TG_CpG_5,col="black",lwd=2,lty=2)
lines(read_position,data$TA_CpG_5,col="black",lwd=2,lty=2)
lines(read_position,data$GA_CpG_5,col="blue",lwd=2,lty=2)
lines(read_position,data$CT_CpG_5,col="red",lwd=2,lty=2)

lines(read_position,data$GA5,col="blue",lwd=2)
lines(read_position,data$CT5,col="red",lwd=2)

plot(read_position,read_position,xlab="Distance from 3' end of sequence read",ylab="Mismatch frequency",cex.axis=1,cex.lab=1.5,cex=3,type="n",ylim=c(0,0.50),xlim=c(max(read_position),0))
lines(read_position,data$CG3,col="black",lwd=2)
lines(read_position,data$CA3,col="black",lwd=2)
lines(read_position,data$GT3,col="black",lwd=2)
lines(read_position,data$GC3,col="black",lwd=2)
lines(read_position,data$AG3,col="black",lwd=2)
lines(read_position,data$AT3,col="black",lwd=2)
lines(read_position,data$AC3,col="black",lwd=2)
lines(read_position,data$TC3,col="black",lwd=2)
lines(read_position,data$TG3,col="black",lwd=2)
lines(read_position,data$TA3,col="black",lwd=2)

lines(read_position,data$CG_CpG_3,col="black",lwd=2,lty=2)
lines(read_position,data$CA_CpG_3,col="black",lwd=2,lty=2)
lines(read_position,data$GT_CpG_3,col="black",lwd=2,lty=2)
lines(read_position,data$GC_CpG_3,col="black",lwd=2,lty=2)
lines(read_position,data$AG_CpG_3,col="black",lwd=2,lty=2)
lines(read_position,data$AT_CpG_3,col="black",lwd=2,lty=2)
lines(read_position,data$AC_CpG_3,col="black",lwd=2,lty=2)
lines(read_position,data$TC_CpG_3,col="black",lwd=2,lty=2)
lines(read_position,data$TG_CpG_3,col="black",lwd=2,lty=2)
lines(read_position,data$TA_CpG_3,col="black",lwd=2,lty=2)
lines(read_position,data$GA_CpG_3,col="blue",lwd=2,lty=2)
lines(read_position,data$CT_CpG_3,col="red",lwd=2,lty=2)

lines(read_position,data$GA3,col="blue",lwd=2)
lines(read_position,data$CT3,col="red",lwd=2)

plot.new()
legend("center",c("C->T","G->A","CpG->TpG","CpG->CpA","All others"),lty=c(1,1,2,2,1),lwd=c(2,2,2,2,2),col=c("red","blue","red","blue","black"),bty="b",cex=2,bg="white")

dev.off()
