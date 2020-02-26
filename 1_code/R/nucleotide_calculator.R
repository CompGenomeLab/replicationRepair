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
        "initiation.zones.repdomains.hela_windows_201_100.bed", sep = ""), 
  header=FALSE)
names(dataframe) <- c("chr", "start", "end", "id", "value", "strand")

id <- dataframe[4]
chr <- dataframe[1]
df_start <- dataframe[2]
df_end <- dataframe[3]

finaltable <- data.frame(id = unique(dataframe$id), patC = 0, patG = 0, 
                         length = 0, strand = "+")

pattern1<-DNAStringSet("C")
pattern2<-DNAStringSet("G")

row <- 1

while (row <= nrow(id)) {
  
  f_row <- (1:nrow(finaltable))[finaltable$id == id[row,]]

  Dnastring <- DNAStringSet(Hsapiens[[as.character(chr[row,])]], df_start[row,], 
                          df_end[row,])
  
  finaltable$patC[f_row] <- finaltable$patC[f_row] + 
    vcountPDict(pattern1, Dnastring)
  finaltable$patG[f_row] <- finaltable$patG[f_row] + 
    vcountPDict(pattern2, Dnastring)
  
  finaltable$length[f_row] <- finaltable$length[f_row] + 
    df_end[row,] - df_start[row,]

  if (row %% 201 == 0){ 
    print(row) 
  }
  
  row <- row + 1 
}

finaltable$freqC <- finaltable$patC / finaltable$length * 100
finaltable$freqG <- finaltable$patG / finaltable$length * 100

backup_just_in_case<-finaltable
setwd("~/Desktop")
write.table(finaltable, "InZones_nuc_content.txt", sep="\t", 
            row.names = FALSE, quote = FALSE)



df <- finaltable


for (rearrange in 1) {  
  
  #### set variables ####
  
  dateout <- Sys.Date()
  dateout <- format(dateout, format = "[%Y.%m.%d]")
  
  window_number <- 201
  if (window_number %% 2 == 0) {
    half_window <- window_number / 2
  } else {
    half_window <- (window_number - 1) / 2 + 1
  }
  
  #### rearrange ####
  
  # window numbers separated from dataset names
  window <- data.frame(str_split_fixed(df$id, "_", -1))
  
  df$windows <- as.numeric(levels(
    window[ , ncol(window)]))[window[ , ncol(window)]] - half_window
  
  df$id <- "Initiation Zones"
  
  df$repdomains <- window[,3]
  
  rm(window, rearrange)
  
  
}


#### plot ####

#### filter the data ####

#d <- filter(df, repdomains == "DTZ") # filter

df$repdomains <- factor(df$repdomains, levels = c("UTZ", "ERD", "DTZ", "LRD"))

#### add plot format #### 

sourcePath <- paste("~/Documents/myprojects/replicationRepair/1_code/R/")

source(paste(sourcePath, "4_plot_format.R", sep = ""))

#### create the plot ####

p <- ggplot() + 
  geom_line(data = df, aes(x = windows, y = freqC, color = "freqC")) + 
  geom_line(data = df, aes(x = windows, y = freqG, color = "freqG")) + 
  facet_grid(~repdomains) +
  xlab(windows_lab) + ylab("Cytosine Frequency (%)") +
  scale_x_continuous(limits = c(-100, 100), 
                     breaks = c(-100, 0, 100), 
                     labels = c("-10 kb", "Initiation Zones", "+10 kb")) + 
  scale_color_manual(values = c(freqC = "#0571b0", 
                                freqG = "#ca0020"), 
                     labels = c("+", "-")) +
  labs(color = "Strands") +
  labs(title = paste("Cytosine Content of Initiation Zones"), 
       subtitle = paste("20 kb sliding windows at 100 bp intervals", sep = ""))

p <- p + p_format  # add plot format

p # visualize

#### save the plot ####

figurePath <- paste("~/Documents/myprojects/replicationRepair/", 
                    "4_output/gitignore/InZones/", sep = "")

figureName <- paste(dateout, "C_Content_of_Initiation_Zones_with_repdomains.pdf", sep = "")

ggsave(path = figurePath, filename = figureName, 
       width = 297, height = 210, units = "mm")

figurePNG <- sub(".pdf", ".png", figureName)

ggsave(path = figurePath, filename = figurePNG, 
       width = 297, height = 210, units = "mm")

#### save the figure info ####

#dataInfo <- paste(date, "final_report_", fr_name, "_info.TXT", sep = "")

figureInfo <- "nucleotide_calculator.R"

source(paste(sourcePath, "4_figure_info.R", sep = ""))







