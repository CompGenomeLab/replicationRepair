#### Packages and Libraries ####

library(argparser)
library(stringr)
library(ggpubr)
library(ggplot2)
library(dplyr)
library(reshape2)
library(patchwork)
set.seed(1) 

######## Arguments ##########
p <- arg_parser("producing the figure 2")
p <- add_argument(p, "--windowed", help="windowed (20kb) region file with read counts")
p <- add_argument(p, "--windowed_sim", help="windowed (20kb) region file with simulated read counts")
p <- add_argument(p, "--noW", help="Replication domains with read counts")
p <- add_argument(p, "--noW_sim", help="Replication domains with simulated read counts")
p <- add_argument(p, "--data_prefix", help="name prefix of the dataframes that generate the plots")
p <- add_argument(p, "-o", help="output")

# Parse the command line arguments
argv <- parse_args(p)

sample_csv_pA <- argv$windowed

sample_csv_pB <- argv$noW

sample_csv_pA_sim <- argv$windowed_sim

sample_csv_pB_sim <- argv$noW_sim

# name of the sample csv file for plot A
#sample_csv_pA <- paste("/Users/azgarian/Desktop/replication_final/",
#                       "final_reports_hg19_repdomains_uv_mean0.5_windows_201_10000.txt", 
#                       sep = "")

#sample_csv_pA_sim <- paste("/Users/azgarian/Desktop/replication_final/",
#                       "final_reports_sim_hg19_repdomains_uv_mean0.5_windows_201_10000.txt", 
#                       sep = "")

# name of the sample csv file for plot B
#sample_csv_pB <- paste("/Users/azgarian/Desktop/replication_final/",
#                       "final_reports_noWindows_hg19_repdomains_hela.txt", 
#                       sep = "")

#sample_csv_pB_sim <- paste("/Users/azgarian/Desktop/replication_final/",
#                           "final_reports_noWindows_sim_hg19_repdomains_hela.txt", 
#                           sep = "")

#### Default Plot Format ####

source("workflow/scripts/plot_format.R")


#### Functions ####

source("workflow/scripts/functions.R")

real_sim_pA <- function ( df, sim_rep = "_" ){
  sim <- filter(df, type == "sim" & replicate == sim_rep)
  real <- filter(df, type == "real")
  real_list <- split(real, real$replicate)
  df_real_sim <- real_list[[1]][0, ]
  
  for ( real_rep in 1:length(real_list) ){
    temp <- real_list[[real_rep]]
    temp <- rbind(temp, sim) 
    temp <- dcast(temp, chromosomes + start_position + end_position + 
                    dataset + score + dataset_strand + product + phase + windows +
                    time_after_exposure + sample_strand ~ type, 
                  value.var = "xr_ds")
    temp <- cbind(temp, real_list[[real_rep]]["replicate"])
    df_real_sim <- rbind(df_real_sim, temp)
  }
  
  df_real_sim$real_sim <- df_real_sim$real / df_real_sim$sim
  df_real_sim <- select(df_real_sim, -c("real", "sim"))
  
  return(df_real_sim)
}

#### Variables ####

rep <- "_" 
if (rep == "A"){ 
  fig_name = "~/Desktop/fig2.png" 
  ann_text_ERD <- data.frame(phase = 0.7, real_sim = 4.8,
                             dataset = factor("ERD",levels = c("ERD","LRD")),
                             product = factor("64_PP", 
                                              levels = c("64_PP", "CPD")))
  lab_erd <- ""
  ann_text_LRD <- data.frame(phase = 0.9, real_sim = 4.8,
                             dataset = factor("LRD",levels = c("ERD","LRD")),
                             product = factor("64_PP", 
                                              levels = c("64_PP", "CPD")))
  lab_lrd <- ""
} else if (rep == "B"){ 
  fig_name = "~/Desktop/supfig4.png" 
  ann_text_ERD <- data.frame(phase = 0.7, real_sim = 4.8,
                             dataset = factor("ERD",levels = c("ERD","LRD")),
                             product = factor("64_PP", 
                                              levels = c("64_PP", "CPD")))
  lab_erd <- ""
  ann_text_LRD <- data.frame(phase = 0.7, real_sim = 4.7,
                             dataset = factor("LRD",levels = c("ERD","LRD")),
                             product = factor("64_PP", 
                                              levels = c("64_PP", "CPD")))
  lab_lrd <- ""
} else if (rep == "_"){ 
  fig_name = "~/Desktop/fig2_comb.png" 
  ann_text_ERD <- data.frame(phase = 0.7, real_sim = 4.7,
                             dataset = factor("ERD",levels = c("ERD","LRD")),
                             product = factor("64_PP", 
                                              levels = c("64_PP", "CPD")))
  lab_erd <- ""
  ann_text_LRD <- data.frame(phase = 0.7, real_sim = 4.7,
                             dataset = factor("LRD",levels = c("ERD","LRD")),
                             product = factor("64_PP", 
                                              levels = c("64_PP", "CPD")))
  lab_lrd <- ""
} 

#### Main ####

# for plot A
pA_sample_df <- read.delim( sample_csv_pA, header = F )

colnames(pA_sample_df) <- c("chromosomes", "start_position", "end_position", 
                            "dataset", "score", "dataset_strand", "counts", 
                            "sample_names", "file_names", "layout", "cell_line", 
                            "product", "method", "uv_exposure", "treatment", 
                            "phase", "time_after_exposure", "replicate", 
                            "project", "sample_source", "sample_strand", 
                            "mapped_reads", "RPKM")

pA_df_rr <- repair_rate( pA_sample_df, ds_rep = "_" )
pA_df_rr_org <- window_numbering( pA_df_rr, 4, 101 )
pA_df_rr_org$dataset <- gsub("_.*", "", pA_df_rr_org$dataset)

pA_sample_df_sim <- read.delim( sample_csv_pA_sim, header = F )

colnames(pA_sample_df_sim) <- c("chromosomes", "start_position", "end_position", 
                            "dataset", "score", "dataset_strand", "counts", 
                            "sample_names", "file_names", "layout", "cell_line", 
                            "product", "method", "uv_exposure", "treatment", 
                            "phase", "time_after_exposure", "replicate", 
                            "project", "sample_source", "sample_strand", 
                            "mapped_reads", "RPKM")

pA_df_rr_sim <- repair_rate( pA_sample_df_sim, ds_rep = "_" )
pA_df_rr_org_sim <- window_numbering( pA_df_rr_sim, 4, 101 )
pA_df_rr_org_sim$dataset <- gsub("_.*", "", pA_df_rr_org_sim$dataset)

pA_df_rr_org$type <- "real"
pA_df_rr_org_sim$type <- "sim"

pA_df_rr_comb <- rbind(pA_df_rr_org, pA_df_rr_org_sim)

pA_df_rr_rs <- real_sim_pA(pA_df_rr_comb)


# for plot B
pB_sample_df <- read.delim( sample_csv_pB, header = F )

colnames(pB_sample_df) <- c("chromosomes", "start_position", "end_position", 
                            "dataset", "score", "dataset_strand", "counts", 
                            "sample_names", "file_names", "layout", "cell_line", 
                            "product", "method", "uv_exposure", "treatment", 
                            "phase", "time_after_exposure", "replicate", 
                            "project", "sample_source", "sample_strand", 
                            "mapped_reads", "RPKM")

pB_sample_df <- rmv_low_counted_regions( pB_sample_df )
pB_df_rr <- repair_rate( pB_sample_df, ds_rep = "_" )

pB_sample_df_sim <- read.delim( sample_csv_pB_sim, header = F )

colnames(pB_sample_df_sim) <- c("chromosomes", "start_position", "end_position", 
                            "dataset", "score", "dataset_strand", "counts", 
                            "sample_names", "file_names", "layout", "cell_line", 
                            "product", "method", "uv_exposure", "treatment", 
                            "phase", "time_after_exposure", "replicate", 
                            "project", "sample_source", "sample_strand", 
                            "mapped_reads", "RPKM")

pB_sample_df_sim <- rmv_low_counted_regions( pB_sample_df_sim )
pB_df_rr_sim <- repair_rate( pB_sample_df_sim, ds_rep = "_" )

pB_df_rr$type <- "real"
pB_df_rr_sim$type <- "sim"

pB_df_rr_comb <- rbind(pB_df_rr, pB_df_rr_sim)

pB_df_rr_rs <- dcast(pB_df_rr_comb, chromosomes + start_position + end_position + 
                       dataset + score + dataset_strand + product + phase +
                       time_after_exposure + sample_strand + replicate ~ type, 
                     value.var = "xr_ds")

pB_df_rr_rs$real_sim <- pB_df_rr_rs$real / pB_df_rr_rs$sim
pB_df_rr_rs <- select(pB_df_rr_rs, -c("real", "sim"))

#### Filtering Samples ####

pA1_data <- filter(pA_df_rr_rs, phase != "async", dataset == "ERD", 
                   replicate == rep, time_after_exposure == "12")

pA2_data <- filter(pA_df_rr_rs, phase != "async", dataset == "LRD", 
                   replicate == rep, time_after_exposure == "12")

pB_data <- filter(pB_df_rr_rs, dataset == "ERD" | dataset == "LRD", 
                 replicate == rep, phase != "async", 
                 time_after_exposure == "12", real_sim != "NA")

pB_data$real_sim <- as.numeric(pB_data$real_sim)

#### Plot A.1 ####

write.table(pA1_data, file = paste0(argv$data_prefix, "A1.csv"), 
            quote = FALSE, row.names = FALSE, sep = ",")

# create the plot
p.A.1 <- ggplot(pA1_data, aes(x = windows, y = log2(real_sim))) + 
  geom_hline(yintercept = 0, linetype="dashed", color="red") +
  geom_line(aes(color = phase, linetype = sample_strand)) +
  facet_grid(~product~time_after_exposure~dataset,
             labeller = labeller(product = product_labs, 
                                 time_after_exposure = taex_labs)) +
  xlab("Relative Position (kb)") + ylab("n. Repair Rate (log2)") +
  scale_x_continuous(limits = c(-101, 101), 
                     breaks = c(-101, 0, 101), 
                     labels = c("-1000", "0", "+1000")) + 
  scale_color_manual(name = "Phase", 
                     label = c("Early S Phase", "Late S Phase", "Asyncronized"), 
                     values = phase_colors) + 
  labs(color = "Phase", linetype = "Strands") +
  ylim(-2, 1) +
  guides(linetype = "none")

# adding and overriding the default plot format
p.A.1 <- p.A.1 + p_format + 
  theme(strip.text.y = element_text(size=0, margin = margin(0, 0, 0, 0)),
        legend.position = "right",
        axis.text.x = element_text(hjust=c(0.1, 0.5, 0.9)))


#### Plot A.2 ####

write.table(pA2_data, file = paste0(argv$data_prefix, "A2.csv"), 
            quote = FALSE, row.names = FALSE, sep = ",")

# create the plot
p.A.2 <- ggplot(pA2_data, aes(x = windows, y = log2(real_sim))) +
  geom_hline(yintercept = 0, linetype="dashed", color="red") +
  geom_line(aes(color = phase, linetype = sample_strand)) +
  facet_grid(~product~time_after_exposure~dataset,
             labeller = labeller(product = product_labs, 
                                 time_after_exposure = taex_labs)) +
  xlab("Relative Position (kb)") + ylab("n. Repair Rate (log2)") +
  scale_x_continuous(limits = c(-101, 101), 
                     breaks = c(-101, 0, 101), 
                     labels = c("-1000", "0", "+1000")) + 
  scale_color_manual(name = "Phase", 
                     label = c("Early S Phase", "Late S Phase", "Asyncronized"), 
                     values = phase_colors) + 
  labs(color = "Phase", linetype = "Strands") +
  ylim(-2, 1) +
  guides(linetype = "none")

# adding and overriding the default plot format
p.A.2 <- p.A.2 + p_format + 
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        strip.text.y = element_text(size=0, margin = margin(0, 0, 0, 0)),
        legend.position = "right",
        axis.text.x = element_text(hjust=c(0.1, 0.5, 0.9)))

#### Plot B ####

write.table(pB_data, file = paste0(argv$data_prefix, "B.csv"), 
            quote = FALSE, row.names = FALSE, sep = ",")

# create the plot
p.B <- ggplot(pB_data, aes(x = phase, y = log2(real_sim))) + 
  geom_boxplot(aes(fill = phase), outlier.shape = NA) +
  geom_hline(yintercept = 0, linetype="dashed", color="red") +
  facet_grid(product~time_after_exposure~dataset, 
             labeller = labeller(product = product_labs, 
                                 time_after_exposure = taex_labs, 
                                 replicate = rep_labs)) + 
  xlab(phase_lab) + ylab("") +
  scale_x_discrete(labels = c("Early S", "Late S")) +
  scale_fill_manual(name="", 
                    values = phase_colors) +
  guides(fill = "none") 

# adding and overriding the default plot format
# p values calculated
p.B <- p.B + p_format + 
  theme(axis.title.y=element_blank()) +
  stat_compare_means( comparisons = list( c("early", "late") ), 
                      label.y = 2) +
  geom_text(data = ann_text_ERD, label = lab_erd) +
  geom_text(data = ann_text_LRD, label = lab_lrd)

#### Combining Plots with Patchwork ####
p_final <- (p.A.1 + labs(title="A")) + p.A.2 + (p.B  + labs(title="B")) +
  plot_layout(guides = "collect") &
  theme(plot.tag = element_text(size = 12, face="bold"),
        legend.position = 'bottom', 
        plot.title = element_text(hjust = -0.2, 
                                  size = 12, face="bold"))

ggsave(argv$o, width = 22, height = 18, units = "cm")
