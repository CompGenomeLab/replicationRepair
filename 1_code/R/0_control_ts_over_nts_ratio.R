#### library ####

library("ggplot2")
library(grid)
library(gridExtra)
library(dplyr)
library(tidyr)
library(reshape2)

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

# trash <- filter(all_ts_nts, TSoverNTS == Inf | is.nan(TSoverNTS))
# 
# all_ts_nts <- all_ts_nts[!(all_ts_nts$gene_id %in% trash$gene_id), ]  

mynames <- c("cell","product","method","treatment","phase","exposure","replicate")

fr <- separate(all_ts_nts, c("name"), into = mynames, sep = "_")

ds_a <- filter(fr, method == "DS" & replicate == "A")

fr_xr <- filter(fr, method == "XR")

fr_list <- split(fr_xr, fr_xr$replicate)

fr_xr_ds <- fr_list[[1]][0, ]

for (xr in 1:length(fr_list)) {
  
  temp <- fr_list[[xr]]
  
  temp <- rbind(temp, ds_a)
  
  temp <- dcast(temp, chromosomes + s_point + e_point + 
                  gene_id + strand + cell + product + phase + 
                  exposure ~ method, 
                value.var = "TSoverNTS")
  
  temp <- cbind(temp, fr_list[[xr]]["replicate"])
  
  fr_xr_ds <- rbind(fr_xr_ds, temp)
}

fr_xr_ds$xr_ds <- fr_xr_ds$XR / fr_xr_ds$DS

fr_xr_ds = select(fr_xr_ds, -c("DS", "XR"))

rm(ds_a, fr_xr, fr_list, xr, temp)

fr_xr_ds$name <- paste(fr_xr_ds$cell, fr_xr_ds$product, fr_xr_ds$phase, 
                       fr_xr_ds$exposure)


#### plot ####

p <- ggplot(fr_xr_ds, aes(x = name, y = log2(xr_ds))) + 
  geom_boxplot(outlier.shape = NA) + 
  xlab("Sample Names") + ylab("log2 normalized TS/NTS of Repair Rate") +
  ylim(-1, 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") 

#### theme ####

p <- p + theme_light() +
  theme(axis.title.x = element_text(size = 16, face = "bold"),
        axis.title.y = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 16, angle = -60, vjust = 0.5, 
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

ggsave('TS_NTS_repair_rate_log2norm.png', width = 397, height = 210, units = "mm") 


