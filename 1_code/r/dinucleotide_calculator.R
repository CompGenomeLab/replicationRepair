library(ggplot2)
library(reshape2)
library(dplyr)
library(plyr)
library(BSgenome)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
library(Biostrings)

#ADEBALI LAB - SABANCI UNIVERSITY
#EXON REPAIR PROJECT
#Ogün Adebali, Cem Azgari, Berk Turhan, Zeynep Kilinç, Defne Çirci

#This R script aims calculate the count of each dinucleotide in each exonic or intronic region according to input file

###Inputs: One file with all exonic regions or intronic regions in bed format
###Outputs: One file with dinucleotide counts of specified regions

#Dinucleotide Calculator Code#

#Initilize the files folder 
setwd("~/Desktop")

#Enter the file to calculate dinucleotide counts(Intron file or Exon file)
dataframe<-read.table(
  paste("~/Documents/myprojects/replicationRepair/0_data/gitignore/InZones/", 
        "initiation.zones.repdomains.hela_windows_201_100.bed", sep = ""), 
  header=FALSE)
names(dataframe) <- c("chr", "start", "end", "id", "value", "strand")
finaltable <- dataframe[c(1:201),c(4,6)]

#This data frame may indicate exon or intron data depending on the entry file
region<-dataframe
region$id<- as.character.factor(region$id)

finaltable$strand <- "+"
finaltable$TT <- 0
finaltable$TC <- 0 
finaltable$CT <- 0
finaltable$CC <- 0

#Calculation of dinucleotide counts
for  (row in 1:(length(dataframe$chr))){
  
  modrow <- 201 %% row
  if (201 %% row == 0){ modrow = 201 } 
  
  chromosom<-as.character(region[row,1])
  Dnastring<-DNAStringSet(Hsapiens[[chromosom]],region$start[row],region$end[row])
  
  #The non transcribed strands reverse complementary sequence is taken so that
  #all dinucleotide counts are calculated according to transcribed strand
  if(as.character(region$strand[row])=="-"){
    Dnastring=reverseComplement(Dnastring)}
  
  pattern1<-DNAStringSet("TT")
  pattern2<-DNAStringSet("TC")
  pattern3<-DNAStringSet("CT")
  pattern4<-DNAStringSet("CC")
  
  finaltable$TT[modrow] <- finaltable$TT[modrow] + vcountPDict(pattern1,Dnastring)
  finaltable$TC[modrow] <- finaltable$TC[modrow] + vcountPDict(pattern2,Dnastring)
  finaltable$CT[modrow] <- finaltable$CT[modrow] + vcountPDict(pattern3,Dnastring)
  finaltable$CC[modrow] <- finaltable$CC[modrow] + vcountPDict(pattern4,Dnastring)
  
  if (row %% 1000000 == 0){ 
    print(row/201)
    }
}

#In case of a problem save the data frame
backup_just_in_case<-finaltable

write.table(finaltable, "InZones_plus_dinuc_content.txt", sep="\t", 
            row.names = FALSE, quote = FALSE)
