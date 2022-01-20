#### Packages and Libraries ####

library(argparser)
library(stringr)
library(ggpubr)
library(ggplot2)
library(dplyr)
library(reshape2)
library(patchwork)

######## Arguments ##########
p <- arg_parser("producing damage and repair signals that are normalized by simulations.")
p <- add_argument(p, "--df", help="region file with read counts")
p <- add_argument(p, "--df_sim", help="region file with simulated read counts")
p <- add_argument(p, "-o", help="output")

# Parse the command line arguments
argv <- parse_args(p)

# name of the sample csv file for plot A
sample_pA <- argv$df

sample_pA_sim <- argv$df_sim

#### Default Plot Format ####

source("workflow/scripts/plot_format.R")


#### Fuctions ####

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
pA_df_org$repDomain <- gsub("iz_", "", pA_df_org$dataset)
pA_df_org$repDomain <- gsub("_.*", "", pA_df_org$repDomain)

pA_sample_df_sim <- read.delim( sample_pA_sim, header = F )
colnames(pA_sample_df_sim) <- c("chromosomes", "start_position", "end_position", 
                                "dataset", "score", "dataset_strand", "counts", 
                                "sample_names", "file_names", "layout", "cell_line", 
                                "product", "method", "uv_exposure", "treatment", 
                                "phase", "time_after_exposure", "replicate", 
                                "project", "sample_source", "sample_strand", 
                                "mapped_reads", "RPKM")

pA_df_org_sim <- window_numbering( pA_sample_df_sim, 4, 101 )
pA_df_org_sim$repDomain <- gsub("iz_", "", pA_df_org_sim$dataset)
pA_df_org_sim$repDomain <- gsub("_.*", "", pA_df_org_sim$repDomain)

colnames(pA_df_org)[23] <- "real"
colnames(pA_df_org_sim)[23] <- "sim"

pA_df_org2 <- pA_df_org[ -c(1:7,9,10,14,15,19,20,22) ] 
pA_df_org2_sim <- pA_df_org_sim[ -c(1:7,9,10,14,15,19,20,22) ] 

pA_df_rs <- merge(pA_df_org2, pA_df_org2_sim, by = c("repDomain", "sample_names", "cell_line", 
                                                   "product", "method", "phase", 
                                                   "time_after_exposure", "replicate", 
                                                   "sample_strand", "windows"))
pA_df_rs$xr_ds <- pA_df_rs$real / pA_df_rs$sim

# filtering samples 
pA1_data <- filter(pA_df_rs, phase != "async", 
                   replicate == "_", time_after_exposure == "12",
                   product == "CPD")

pA1_data$time_after_exposure[pA1_data$method == "Damage_seq"] <- "0" 

pA1_data$sample_strand <- factor(
  pA1_data$sample_strand, levels = c("+","-"))

#### Plot A.1 ####

# create the plot
p.A.1 <- ggplot(pA1_data, aes(x = windows, y = log2(xr_ds))) + 
  geom_line(aes(color = sample_strand)) +
  facet_grid(~product~method~time_after_exposure~phase~repDomain,
             labeller = labeller(product = product_labs, 
                                 time_after_exposure = taex_labs,
                                 method = method_labs,
                                 phase = phase_labs)) +
  xlab("Relative Position (kb)") + ylab("Real/Simulation") +
  scale_x_continuous(limits = c(-101, 101), 
                     breaks = c(-101, 0, 101), 
                     labels = c("-10", "0", "+10")) + 
  scale_color_manual(values = strand_colors) +
  labs(color = "Strands") 

# adding and overriding the default plot format
p.A.1 <- p.A.1 + p_format +
  theme(panel.border = element_rect(fill = NA)) 

p.A.1
#### Combining Plots with Patchwork ####


ggsave(argv$o, width = 22, height = 18, units = "cm")

