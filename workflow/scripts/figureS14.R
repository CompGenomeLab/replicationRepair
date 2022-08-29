#### Packages and Libraries ####

library(stringr)
library(ggpubr)
library(ggplot2)
library(dplyr)
library(reshape2)
library(patchwork)
library(grid)
library(argparser)
set.seed(1) 

######## Arguments ##########
p <- arg_parser("producing figure supplementary 14")
p <- add_argument(p, "--df", help="region file with read counts")
p <- add_argument(p, "--df_sim", help="region file with simulated read counts")
p <- add_argument(p, "--data_prefix", help="name prefix of the dataframes that generate the plots")
p <- add_argument(p, "-o", help="output")

# Parse the command line arguments
argv <- parse_args(p)

sample_csv <- argv$df

sample_sim_csv <- argv$df_sim


#### Default Plot Format ####

source("workflow/scripts/plot_format.R")


#### Fuctions ####

source("workflow/scripts/functions.R")


#### Main ####

sample_df <- read.delim( sample_csv, header = F )

colnames(sample_df) <- c("chromosomes", "start_position", "end_position", 
                         "dataset", "score", "dataset_strand", "counts", 
                         "sample_names", "file_names", "layout", "cell_line", 
                         "product", "method", "uv_exposure", "treatment", 
                         "phase", "time_after_exposure", "replicate", 
                         "project", "sample_source", "sample_strand", 
                         "mapped_reads", "RPKM")


sample_df <- window_numbering( sample_df, 4, 101 )
sample_df <- domain_name( sample_df, 1 )
sample_df$dataset <- "Initiation Zones"

sample_df$sample_strand <- factor(
  sample_df$sample_strand, levels = c("+","-"))

sample_df_filt <- filter(sample_df, method == "Damage_seq", phase != "async",
                         product == "CPD")

sample_df_filt_dcast <- 
  dcast(sample_df_filt, dataset + chromosomes + start_position + end_position + 
        score + dataset_strand + product + phase + method + sample_strand + 
        windows + repdomains
        ~ time_after_exposure, value.var = "RPKM")

sample_df_filt_dcast$diff <- 
  (sample_df_filt_dcast$`12` - sample_df_filt_dcast$`120`) / sample_df_filt_dcast$`12`

write.table(sample_df_filt_dcast, file = paste0(argv$data_prefix, ".csv"), 
            quote = FALSE, row.names = FALSE, sep = ",")

p.A <- ggplot(sample_df_filt_dcast, aes(x = windows, y = diff)) +
  geom_vline(xintercept = 0, color = "gray", linetype = "dashed") + 
  geom_line(aes(color = sample_strand)) + 
  facet_grid(~product~phase~repdomains,
             labeller = labeller(product = product_labs, 
                                 phase = phase_labs)) +
  xlab("Relative Position Around Initiation Zones (kb)") + 
  ylab("Normalized Accumulated Repair") +
  scale_x_continuous(limits = c(-101, 101), 
                     breaks = c(-101, 0, 101), 
                     labels = c("-10", "0", "+10")) + 
  scale_color_manual(values = strand_colors) + 
  labs(color = "Strands")

p.A <- p.A + p_format

ggsave( argv$o, width = 22, height = 18, units = "cm" )
