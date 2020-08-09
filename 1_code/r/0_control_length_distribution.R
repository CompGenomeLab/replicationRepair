#### library ####

library("ggplot2")
library(grid)
library(gridExtra)

#### set variables ####

options(scipen=999)

setwd(paste("~/Documents/myprojects/replicationRepair/3_output/gitignore/", 
            "2_sample_control/length_distribution", sep = ""))

only_xr = "_cutadapt_length_distribution.txt"
all_together = "length_distribution.txt"

temp <- list.files(pattern = all_together)

sample_info <- read.csv(paste("~/Documents/myprojects/replicationRepair/", 
                              "0_data/gitignore/", 
                              "project_repair_replication_all_samples.csv", 
                              sep = ""))

#### plot ####

pdf("xr_seq_length_distributions.png", width = 50, height = 30)

p <- list()

for(counter in 1:length(temp)) {
  
  d <- read.delim(temp[[counter]], header = FALSE)
  
  file_name_temp <- gsub("_cutadapt_.*", "", temp[[counter]])
  
  file_name <- gsub("-", "_", file_name_temp)
  
  sample_name <- sample_info$sample_name[grep(file_name, sample_info$file_name)]
  
  colnames(d) <- c("oligomer_length", "counts")
  
  p[[counter]] <- ggplot(d, aes(x = oligomer_length, y = counts, 
                                fill = counts)) + 
    geom_bar(stat = "identity") +
    xlim(15, 35) +
    xlab("Oligomer Length") + ylab("Counts") +
    labs(title="Length Distribution of the Reads", 
         subtitle="Bar Chart",
         caption = sample_name)
  
  
  #### theme ####
  
  p[[counter]] <- p[[counter]] + theme_light() +
    theme(axis.title.x = element_text(size = 16, face = "bold"),
          axis.title.y = element_text(size = 16, face = "bold"),
          axis.text.x = element_text(size = 14, vjust = 0.6),
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
  
  ggsave(gsub(".txt",".pdf", temp[[counter]]), width = 10, height = 8)
}

#### plot save all together ####

do.call(grid.arrange, p)
dev.off()

rm(d, p, sample_info, all_together, counter, file_name, only_xr, sample_name, 
   temp, file_name_temp)
