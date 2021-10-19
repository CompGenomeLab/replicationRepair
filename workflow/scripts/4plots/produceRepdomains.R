library(preprocessCore)

setwd("~/Desktop")
merge<-read.table("merge_RT.txt", header=FALSE)
colnames(merge)<-c("chr","start","end","rep1_T.bg","rep2_T.bg")

merge$rep1_T.bg <- as.double(merge$rep1_T.bg)
merge$rep2_T.bg <- as.double(merge$rep2_T.bg)

merge$avg_T.bg <- (merge$rep1_T.bg + merge$rep2_T.bg) / 2 

merge_values<-as.matrix(merge[,4:ncol(merge)])

ad<-stack(merge[,4:ncol(merge)])$values

norm_data<-normalize.quantiles.use.target(merge_values,ad)
merge_norm<-data.frame(merge[,1:3],norm_data)
colnames(merge_norm)<-colnames(merge)

for(i in 4:ncol(merge_norm))
{write.table(merge_norm[complete.cases(merge_norm[,i]), c(1,2,3,i)],
             gsub(".bg", "_qnorm.bedGraph", colnames(merge_norm)[i]),
             sep="\t",row.names=FALSE, quote=FALSE, col.names=FALSE)}

chrs=grep(levels(as.factor(merge_norm$chr)),pattern="[_YM]",invert=TRUE,value=TRUE)

AllLoess=list()

for(i in 1:(ncol(merge_norm)-3)){
  AllLoess[[i]]=data.frame();
  cat("Current dataset:", colnames(merge_norm)[i+3], "\n");
  for(Chr in chrs){
    RTb=subset(merge_norm, merge_norm$chr==Chr);
    lspan=300000/(max(RTb$start)-min(RTb$start));
    cat("Current chrom:", Chr, "\n");
    RTla=loess(RTb[,i+3] ~ RTb$start, span=lspan);
    RTl=data.frame(c(rep(Chr,times=RTla$n)), RTla$x,
                   merge_norm[which(merge_norm$chr==Chr & merge_norm$start %in%
                                      RTla$x),3],RTla$fitted);
    colnames(RTl)=c("chr","start","end",colnames(RTb)[i+3]);
    if(length(AllLoess[[i]])!=0){
      AllLoess[[i]]=rbind(AllLoess[[i]],RTl)};
    if(length(AllLoess[[i]])==0){
      AllLoess[[i]]=RTl}}}

for(i in 1:length(AllLoess)){write.table(AllLoess[[i]]
  [complete.cases(AllLoess[[i]]),], gsub(".bg","_Loess.bedGraph",
  colnames(AllLoess[[i]]))[4], sep="\t", row.names=FALSE, quote=FALSE,
                                         col.names=FALSE)}


LS <- merge(AllLoess[[1]], AllLoess[[2]], by = c("chr", "start","end"))  
LS <- merge(LS, AllLoess[[3]], by = c("chr", "start","end")) 
cor(LS[,c(4:6)])

RTc = subset(merge, chr =="chr1")    # Subset of raw timing data in chr1
LS1c = subset(AllLoess[[1]], chr =="chr1")    # Subset of smoothed data in chr1
LS2c = subset(AllLoess[[2]], chr =="chr1")
LSAc = subset(AllLoess[[3]], chr =="chr1")

#par(mar=c(2.2,5.1,1,1), mfrow=c(3,1), col="grey", pch=19, cex=0.5, cex.lab=1.8, xaxs="i")

#plot(RTc$rep1_T.bg~RTc$end, ylab="rep1", xaxt="n") # Plot raw data points
#lines(LS1c$rep1_T.bg~LS1c$end, col="blue3", lwd=3) # Overlay loess line

#plot(RTc$rep2_T.bg~RTc$end, ylab="rep2", xaxt="n")
#lines(LS2c$rep2_T.bg~LS2c$end, col="blue3", lwd=3)

#plot(RTc$avg_T.bg~RTc$end, ylab="avg", xaxt="n")
#lines(LSAc$avg_T.bg~LSAc$end, col="blue3", lwd=3)


#### Segmentation

library(DNAcopy)

avg =CNA(merge_norm$avg_T.bg, merge_norm$chr, merge_norm$end, 
            data.type="logratio", sampleid ="untreated_avg")
Seg.avg =segment(avg, nperm=10000, alpha=1e-15, undo.splits="sdundo", 
                 undo.SD=1.5, verbose=2)

#par(ask=T,mar=c(3.1,4.1,1,1))  # Set figure margins; ask before replotting
#plot(Seg.avg, plot.type="c")  # Plot each chromosome separately
#plot(Seg.avg, plot.type="s")  # Plot overview of all chromosomes

#plot(subset(Seg.avg,chromlist="chr2"), pch=19, pt.cols=c("gray","gray"),
#     xmaploc=T, ylim=c(-3,3)) # Plot a single chromosome with alternate format

#write.table(Seg.avg$output, row.names=F, quote=F, sep="\t")


myAvg =Seg.avg$output   # Extract domain information
myAvg$size =myAvg$loc.end - myAvg$loc.start  # Calculate domain sizes
myAvgEarly =subset(myAvg, myAvg$seg.mean > 0.5) # Create subset of early domains
myAvgLate =subset(myAvg, myAvg$seg.mean < -0.5) # Create subset of late domains
boxplot(myAvgEarly$size, myAvgLate$size, outline=FALSE)  # Distribution of early/late domain sizes

myAvgEarly$repdomains <- "ERD"
myAvgLate$repdomains <- "LRD"
myAvg <- rbind(myAvgEarly, myAvgLate)

myAvg_org <- myAvg[,c(2,3,4,8,6)]
myAvg_org$strand <- "."

write.table(myAvg_org, "repdomains_mean0.5.bed", col.names=F, row.names=F, quote=F, sep="\t")


