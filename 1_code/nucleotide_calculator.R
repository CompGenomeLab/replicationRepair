library(ggplot2)
library(reshape2)
library(dplyr)
library(stringr)
library(plyr)
library(BSgenome)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
library(Biostrings)

dataframe<-read.table(
  paste("~/Documents/myprojects/replicationRepair/0_data/gitignore/InZones/", 
        "initiation.zones.repdomains.hela_windows_201_100.bed", 
        sep = ""), 
  header=FALSE)
names(dataframe) <- c("chr", "start", "end", "id", "value", "strand")

id <- dataframe[4]
chr <- dataframe[1]
df_start <- dataframe[2]
df_end <- dataframe[3]

finaltable <- data.frame(id = unique(dataframe$id), patC = 0, patG = 0, 
                         length = 0, strand = "+")

pattern1<-DNAStringSet("TC")
pattern1.2<-DNAStringSet("CC")
pattern2<-DNAStringSet("GA")
pattern2.2<-DNAStringSet("GG")

row <- 1
DNAStringSet(Hsapiens)
while (row <= nrow(id)) {
  
  f_row <- (1:nrow(finaltable))[finaltable$id == id[row,]]

  Dnastring <- DNAStringSet(Hsapiens[[as.character(chr[row,])]], df_start[row,], 
                          df_end[row,])
  
  finaltable$patC[f_row] <- finaltable$patC[f_row] + 
    vcountPDict(pattern1, Dnastring) + vcountPDict(pattern1.2, Dnastring)
  finaltable$patG[f_row] <- finaltable$patG[f_row] + 
    vcountPDict(pattern2, Dnastring) + vcountPDict(pattern2.2, Dnastring)
  
  finaltable$length[f_row] <- finaltable$length[f_row] + 
    df_end[row,] - df_start[row,]

  if (row %% 1000000 == 0){ 
    print(row) 
  }
  
  row <- row + 1 
}

finaltable$freqC <- finaltable$patC / finaltable$length * 100
finaltable$freqG <- finaltable$patG / finaltable$length * 100

backup_just_in_case<-finaltable
setwd("~/Desktop")
write.table(finaltable, "[2020.12.02]InZones_TC_nuc_content.txt", sep="\t", 
            row.names = FALSE, quote = FALSE)
