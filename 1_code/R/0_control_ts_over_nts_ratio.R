#### library ####

library("ggplot2")
library(grid)
library(gridExtra)

#### set variables ####

options(scipen=999)

setwd(paste("~/Documents/myprojects/replicationRepair/4_output/gitignore/", 
            "2_sample_control/TS_NTS_ratio", sep = ""))

all_together = "TSoverNTScount.txt"

temp <- list.files(pattern = all_together)

sample_info <- read.csv(paste("~/Documents/myprojects/replicationRepair/", 
                              "0_data/gitignore/", 
                              "project_repair_replication_all_samples.csv", 
                              sep = ""))

#### rearrange ####

all_ts_nts <- data.frame()

for(counter in 1:length(temp)) {
  
  d <- read.delim(temp[[counter]], header = FALSE)
  
  file_name_temp <- gsub("_cutadapt_.*", "", temp[[counter]])
  
  file_name <- gsub("-", "_", file_name_temp)
  
  sample_name <- sample_info$sample_name[grep(file_name, sample_info$file_name)]
  
  colnames(d) <- c("chromosomes", "s_point", "e_point", "gene_id", "strand", 
                   "TS", "NTS") 
  
  d$TSoverNTS <- d$TS / d$NTS
  
  d$name <- sample_name
  
  all_ts_nts <- rbind(all_ts_nts, d)
}

trash <- filter(all_ts_nts, TSoverNTS == Inf | TSoverNTS == NaN)

all_ts_nts <- all_ts_nts[!(all_ts_nts$gene_id %in% trash$gene_id), ]  


#### plot ####

p <- ggplot(all_ts_nts, aes(x = name, y = log2(TSoverNTS))) + 
  geom_boxplot(outlier.shape = NA) + 
  xlab("Sample Names") + ylab("log2 normalized TS/NTS Values") +
  ylim(-1, 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") 

#### theme ####

p <- p + theme_light() +
  theme(axis.title.x = element_text(size = 16, face = "bold"),
        axis.title.y = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 16, angle = -90, vjust = 0.5, 
                                   hjust = 0),
        axis.text.y = element_text(size = 14, vjust = 0.1),
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14, angle = 360),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 14),
        plot.caption = element_text(size = 14, 
                                    face = "italic"),
        plot.title = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 16, face = "bold"))

ggsave('TS_NTS_all_log2norm.png', width = 397, height = 210, units = "mm") 


