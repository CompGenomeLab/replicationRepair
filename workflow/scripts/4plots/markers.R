#### Packages and Libraries ####

library(stringr)
library(ggpubr)
library(ggplot2)
library(dplyr)
library(reshape2)
library(patchwork)
library(ggthemes)
library(grid)


#### Variables ####

# name of the sample csv file 
sample_csv <- paste("/Users/azgarian/Desktop/replication_final/", 
                    "final_reports_markers_iz_repdomains_uv_m0.5_hela_windows_201_100_intergenic.txt",
                    sep = "")


#### Default Plot Format ####

source("/Users/azgarian/Documents/myprojects/replicationRepair/workflow/scripts/4plots/4_plot_format.R")


#### Fuctions ####

source("/Users/azgarian/Documents/myprojects/replicationRepair/workflow/scripts/4plots/4_functions.R")


#### Main ####

sample_df <- read.table( sample_csv )

colnames(sample_df) <- c("chromosomes", "start_position", "end_position", "dataset", 
                  "score", "dataset_strand", "counts", "sample_names", 
                  "sample_strand", 
                  "mapped_reads", "RPKM")

df_rr_org <- window_numbering( sample_df, 4, 101 )
df_rr_org <- domain_name( df_rr_org, 1 )
df_rr_org$dataset <- "Initiation Zones"

df_rr_org$sample_strand <- factor(
  df_rr_org$sample_strand, levels = c("+","-"))


# create the plot 
p <- ggplot(df_rr_org, aes(x = windows, y = RPKM, color=sample_names)) +
  geom_vline(xintercept = 0, color = "gray", linetype = "dashed") + 
  geom_line() + 
  facet_grid(~repdomains) +
  xlab("Relative Repair around Initiation Zones (kb)") + ylab("RPKM") +
  scale_x_continuous(limits = c(-100, 100), 
                     breaks = c(-100, 0, 100), 
                     labels = c("-10", "0", "+10"))
p

ggsave("~/Desktop/iz_markers_20kb.pdf", width = 22, height = 18, units = "cm")
