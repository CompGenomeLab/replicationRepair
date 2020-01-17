#### library ####

library(dplyr)
library("ggplot2")

#### set variables ####

dateout <- Sys.Date()
dateout <- format(dateout, format = "[%Y.%m.%d]")

options(scipen = 999)

setwd(paste("~/Documents/My_Projects/Project_Repair_Replication/",
                "Results/1_TextforPlotting/", sep = ""))

#### import and rearrange data ####

data = read.table("[2019.11.14]final_chromosome_damage.txt")

colnames(data) <- c("chromosomes", "start_position", "end_position", "dataset", 
                  "score", "dataset_strand", "counts", "sample_names", 
                  "file_names", "layout", "cell_line", "product", "method", 
                  "uv_exposure", "treatment", "phase", "time_after_exposure", 
                  "replicate", "project", "sample_source", "sample_strand", 
                  "mapped_reads", "RPKM")

data_new = filter(data, chromosomes != "chrM" & chromosomes != "chrY") 
data_new$chromosomes <- factor(data_new$chromosomes, 
                               levels = c("chr1", "chr2", "chr3", "chr4", 
                                          "chr5", "chr6", "chr7", "chr8", 
                                          "chr9", "chr10", "chr11", 
                                          "chr12", "chr13", "chr14", 
                                          "chr15", "chr16", "chr17", 
                                          "chr18", "chr19", "chr20", 
                                          "chr21", "chr22", "chrX"))

#### plot ####

p <- ggplot(data_new,aes(y = RPKM, x = chromosomes)) + 
      geom_boxplot() + 
      facet_grid(~ product ~ replicate ~ time_after_exposure ~ sample_strand 
                 ~ phase) +
      ylab("Normalized Damage Values") + xlab("Chromosomes")

#### theme ####

p <- p +  theme_light() +
          theme(axis.title.x = element_text(size=14, face="bold"),
                axis.title.y = element_text(size=14, face="bold"),
                axis.text.x = element_text(size=12, angle = 60, vjust = 0.5),
                axis.text.y = element_text(size=12,vjust=0.1),
                strip.text = element_text(size = 12),
                legend.title = element_text(size=10, face="bold"),
                legend.text = element_text(size=10))

setwd(paste("~/Documents/My_Projects/Project_Repair_Replication/Results/", 
            "2_sample_control", sep = ""))

ggsave(paste(dateout, 
             'RPKM_Damage_Values_of_all_Samples_at_Every_Chromosome.pdf', 
             sep = ""), 
       width = 17, height = 11)

rm(data, data_new, p, dateout)
