#### Packages and Libraries ####

library(stringr)
library(ggpubr)
library(ggplot2)
library(dplyr)
library(reshape2)
library(patchwork)
library(ggthemes)
library(grid)
library(argparser)

######## Arguments ##########
p <- arg_parser("producing the marker figure")
p <- add_argument(p, "-i", help="input")
p <- add_argument(p, "-o", help="output")

# Parse the command line arguments
argv <- parse_args(p)

#### Variables ####

# name of the sample csv file 
sample_csv <- argv$i
#sample_csv <- paste("/Users/azgarian/Desktop/replication_final/", 
#                    "final_reports_markers_iz_repdomains_uv_m0.5_hela_windows_201_100_intergenic.txt",
#                    sep = "")

#### Default Plot Format ####

source("workflow/scripts/plot_format.R")


#### Fuctions ####

source("workflow/scripts/functions.R")


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
  xlab("Relative Position around Initiation Zones (kb)") + ylab("RPKM") +
  scale_x_continuous(limits = c(-100, 100), 
                     breaks = c(-100, 0, 100), 
                     labels = c("-10", "0", "+10"))
p <- p + p_format

ggsave(argv$o, width = 22, height = 18, units = "cm")
