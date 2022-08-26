#### Packages and Libraries ####

library(argparser)
library(stringr)
library(ggpubr)
library(ggplot2)
library(dplyr)
library(reshape2)
library(patchwork)

######## Arguments ##########
p <- arg_parser("producing the supplementary figure 4")
p <- add_argument(p, "--df", help="region file with read counts")
p <- add_argument(p, "--df_sim", help="region file with simulated read counts")
p <- add_argument(p, "--data_prefix", help="name prefix of the dataframes that generate the plots")
p <- add_argument(p, "-o", help="output")

# Parse the command line arguments
argv <- parse_args(p)

sample <- argv$df

sample_sim <- argv$df_sim

#### Default Plot Format ####

source("workflow/scripts/plot_format.R")


#### Functions ####

source("workflow/scripts/functions.R")


#### Main ####

sample_df <- read.delim( sample, header = F )
colnames(sample_df) <- c("chromosomes", "start_position", "end_position", 
                         "dataset", "score", "dataset_strand", "counts", 
                         "sample_names", "file_names", "layout", "cell_line", 
                         "product", "method", "uv_exposure", "treatment", 
                         "phase", "time_after_exposure", "replicate", 
                         "project", "sample_source", "sample_strand", 
                         "mapped_reads", "RPKM")

df_org <- window_numbering( sample_df, 4, 101 )
df_org$dataset <- gsub("_.*", "", df_org$dataset)

sample_df_sim <- read.delim( sample_sim, header = F )
colnames(sample_df_sim) <- c("chromosomes", "start_position", "end_position", 
                            "dataset", "score", "dataset_strand", "counts", 
                            "sample_names", "file_names", "layout", "cell_line", 
                            "product", "method", "uv_exposure", "treatment", 
                            "phase", "time_after_exposure", "replicate", 
                            "project", "sample_source", "sample_strand", 
                            "mapped_reads", "RPKM")

df_org_sim <- window_numbering( sample_df_sim, 4, 101 )
df_org_sim$dataset <- gsub("_.*", "", df_org_sim$dataset)

colnames(df_org)[23] <- "real"
colnames(df_org_sim)[23] <- "sim"

df_org <- df_org[ -c(2:7,9,10,14,15,19,20,22) ] 
df_org_sim <- df_org_sim[ -c(2:7,9,10,14,15,19,20,22) ] 

df_rs <- merge(df_org, df_org_sim, by = c("dataset", "sample_names", "cell_line", 
                                                   "product", "method", "phase", 
                                                   "time_after_exposure", "replicate", 
                                                   "sample_strand", "windows"))
df_rs$xr_ds <- df_rs$real / df_rs$sim

# filtering samples 
p_data <- filter(df_rs, phase != "async", 
                   replicate == "_", time_after_exposure == "12",
                   product == "CPD")

p_data$time_after_exposure[p_data$method == "Damage_seq"] <- "0" 

#### Plot ####

write.table(p_data, file = paste0(argv$data_prefix, ".csv"), 
            quote = FALSE, row.names = FALSE, sep = ",")

# create the plot
p <- ggplot(p_data, aes(x = windows, y = log2(xr_ds))) + 
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
p <- p + p_format 

ggsave(argv$o, width = 22, height = 18, units = "cm")

