#### library ####

library("ggplot2")
library(dplyr)
library("tidyr")
library(grid)
library(gridExtra)
library(patchwork)

#### set variables ####

options(scipen=999)

setwd("/home/azgarian/Desktop/dinuc")

temp <- list.files(pattern = "_dinucleotideTable.txt")


sample_info <- read.csv(paste("~/Documents/myprojects/replicationRepair/", 
                              "0_data/gitignore/", 
                              "project_repair_replication_all_samples.csv", 
                              sep = ""))

#### rearrange and plot ####


for(counter in 1:length(temp)) {
  
  dinucleotide_table <- read.table(temp[[counter]], header = TRUE)
  
  # rename columns and get their order
  x_order <- c()
  for (i in 2:ncol(dinucleotide_table)) {
    
    colnames(dinucleotide_table)[i] <- c(paste(i - 1, "-", i, sep = ""))
    
    x_order <- c(x_order, paste(i - 1, "-", i, sep = "")) 
    
  }
  
  # add sample names
  file_name_temp <- gsub("_cutadapt_.*", "", temp[[counter]])
  
  file_name <- gsub("-", "_", file_name_temp)
  
  sample_name <- sample_info$sample_name[grep(file_name, sample_info$file_name)]
  
  dt_organized <- dinucleotide_table %>% gather(Position, count, 
                                                2:ncol(dinucleotide_table))
  
  # reorganize
  colnames(dt_organized) <- c("dinucleotides", "positions", "counts")
  
  dt_organized$freq = 100*dt_organized$counts/sum(dinucleotide_table[2])
  
  dt_organized <- filter(dt_organized, dinucleotides == "CC" |
                           dinucleotides == "CT" | dinucleotides == "TC" |
                           dinucleotides == "TT")
  
  dt_organized$dinucleotides = factor(dt_organized$dinucleotides, 
                                      levels = c("CC", "CT", "TC", "TT"))
  
  dt_organized$positions = factor(dt_organized$positions, 
                                  levels = x_order)
  
  # plot
  p <- ggplot(dt_organized, aes(x = positions, y = freq, 
                                           fill = dinucleotides)) + 
    geom_bar(stat = "identity") +
    ylim(0,101) +
    scale_fill_manual(values = c("forestgreen", "royalblue3",
                                 "gold1", "mediumvioletred")) +
    ylab("Frequency (%)") + xlab("Position in oligomers") +
    labs(title = "Dinucleotide Content of Oligomers", 
         subtitle = "Categorywise Bar Chart", 
         caption = temp[[counter]])
  
  #### theme ####
  
  p <- p + theme_light() +
    theme(axis.title.x = element_text(size = 16, face = "bold"),
          axis.title.y = element_text(size = 16, face = "bold"),
          axis.text.x = element_text(size = 14, vjust = 0.6, angle = 65),
          axis.text.y = element_text(size = 14, vjust = 0.1),
          strip.text.x = element_text(size = 14),
          strip.text.y = element_text(size = 14, angle = 360),
          legend.title = element_text(size = 14, face = "bold"),
          legend.text = element_text(size = 14),
          plot.caption = element_text(size = 14, 
                                      face = "italic"),
          plot.title = element_text(size = 16, face = "bold"),
          plot.subtitle = element_text(size = 16, face = "bold"))
  
  #### plot save 1 by 1 ####
  
  ggsave(gsub(".txt",".pdf", temp[[counter]]), width = 10, 
         height = 8)
  
  if (counter %% 2 == 1) {
    p_total <- p + plot_layout(guides = "collect", tag_level = 'new', ncol = 2) +  
      plot_annotation(tag_levels = c('A'), tag_suffix = ':')
  } else {
    p_total <- p_total + p 
    ggsave(paste("compare_real_to_sim_", 
                 gsub(".txt",".png", temp[[counter]]), sep = ""), 
           width = 497, height = 210, units = "mm")
  }
}




