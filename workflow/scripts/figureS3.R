#### Packages and Libraries ####

library(stringr)
library(ggpubr)
library(ggplot2)
library(dplyr)
library(reshape2)
library(patchwork)
library(ggthemes)

######## Arguments ##########
p <- arg_parser("producing the supplementary figure 3")
p <- add_argument(p, "--df", help="region file with read counts")
p <- add_argument(p, "--df_sim", help="region file with simulated read counts")
p <- add_argument(p, "-o", help="output")

# Parse the command line arguments
argv <- parse_args(p)

# name of the sample csv file for plot A
#sample_pA <- paste("/Users/azgarian/Desktop/replication_final/",
#                   "final_reports_hg19_repdomains_uv_mean0.5_windows_201_10000.txt", 
#                   sep = "")

#sample_pA_sim <- paste("/Users/azgarian/Desktop/replication_final/",
#                   "final_reports_sim_hg19_repdomains_uv_mean0.5_windows_201_10000.txt", 
#                   sep = "")

sample_pA <- argv$df

sample_pA_sim <- argv$df_sim

#### Default Plot Format ####

source("workflow/scripts/plot_format.R")


#### Functions ####

source("workflow/scripts/functions.R")


#### Main ####

# for plot A
pA_sample_df <- read.delim( sample_pA, header = F )
colnames(pA_sample_df) <- c("chromosomes", "start_position", "end_position", 
                         "dataset", "score", "dataset_strand", "counts", 
                         "sample_names", "file_names", "layout", "cell_line", 
                         "product", "method", "uv_exposure", "treatment", 
                         "phase", "time_after_exposure", "replicate", 
                         "project", "sample_source", "sample_strand", 
                         "mapped_reads", "RPKM")

pA_df_org <- window_numbering( pA_sample_df, 4, 101 )
pA_df_org$dataset <- gsub("_.*", "", pA_df_org$dataset)

pA_sample_df_sim <- read.delim( sample_pA_sim, header = F )
colnames(pA_sample_df_sim) <- c("chromosomes", "start_position", "end_position", 
                            "dataset", "score", "dataset_strand", "counts", 
                            "sample_names", "file_names", "layout", "cell_line", 
                            "product", "method", "uv_exposure", "treatment", 
                            "phase", "time_after_exposure", "replicate", 
                            "project", "sample_source", "sample_strand", 
                            "mapped_reads", "RPKM")

pA_df_org_sim <- window_numbering( pA_sample_df_sim, 4, 101 )
pA_df_org_sim$dataset <- gsub("_.*", "", pA_df_org_sim$dataset)

colnames(pA_df_org)[23] <- "real"
colnames(pA_df_org_sim)[23] <- "sim"

pA_df_org <- pA_df_org[ -c(2:7,9,10,14,15,19,20,22) ] 
pA_df_org_sim <- pA_df_org_sim[ -c(2:7,9,10,14,15,19,20,22) ] 

pA_df_rs <- merge(pA_df_org, pA_df_org_sim, by = c("dataset", "sample_names", "cell_line", 
                                                   "product", "method", "phase", 
                                                   "time_after_exposure", "replicate", 
                                                   "sample_strand", "windows"))
pA_df_rs$xr_ds <- pA_df_rs$real / pA_df_rs$sim

# filtering samples 
pA1_data <- filter(pA_df_rs, phase != "async", 
                   replicate == "_", time_after_exposure == "12",
                   product == "CPD")

#### Plot A.1 ####

# create the plot
p.A.1 <- ggplot(pA1_data, aes(x = windows, y = log2(xr_ds))) + 
  geom_smooth(aes(color = phase, linetype = sample_strand), se=FALSE) +
  facet_grid(~product~time_after_exposure~method~dataset,
             labeller = labeller(product = product_labs, 
                                 time_after_exposure = taex_labs,
                                 method = method_labs)) +
  xlab("Relative Position (kb)") + ylab("Real/Simulation") +
  scale_x_continuous(limits = c(-101, 101), 
                     breaks = c(-101, 0, 101), 
                     labels = c("-1000", "0", "+1000")) + 
  scale_color_manual(name = "Phase", 
                     label = c("Early Phase", "Late Phase", "Asyncronized"), 
                     values = phase_colors) + 
  labs(color = "Phase", linetype = "Strands") 

# adding and overriding the default plot format
p.A.1 <- p.A.1 + p_format 

p.A.1
#### Combining Plots with Patchwork ####


ggsave(argv$o, width = 22, height = 18, units = "cm")

