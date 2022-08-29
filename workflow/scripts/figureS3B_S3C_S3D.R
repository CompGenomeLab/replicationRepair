#### Packages and Libraries ####

library(argparser)
library(futile.logger)
library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(patchwork)
library("viridis")
set.seed(1) 

######## Arguments ##########
p <- arg_parser("producing the figure S3B-D")
p <- add_argument(p, "--dinuc_ds_cpd_120", help="dinuc content of a ds sample")
p <- add_argument(p, "--dinuc_xr_cpd_120", help="dinuc content of a xr sample")
p <- add_argument(p, "--len_xr_cpd_120", help="length dist of a xr sample")
p <- add_argument(p, "--tss", help="input tss file")
p <- add_argument(p, "--tes", help="input tes file")
p <- add_argument(p, "--data_prefix", help="name prefix of the dataframes that generate the plots")
p <- add_argument(p, "-o", help="output")
p <- add_argument(p, "--log", help="log file")

# Parse the command line arguments
argv <- parse_args(p)

# log file
flog.appender(appender.file(argv$log))

#### Default Plot Format ####

source("workflow/scripts/plot_format.R")

#### Fuctions ####
flog.info("Importing the functions...") 

source("workflow/scripts/functions.R")

dinuc <- function( dinucleotide_table, method ){
  # rename columns and get their order
  x_order <- c()
  
  if ( method == "XR" ){
    for (i in 2:ncol(dinucleotide_table)) {
      
      colnames(dinucleotide_table)[i] <- c(paste(i - 1, "-", i, sep = ""))
      
      x_order <- c(x_order, paste(i - 1, "-", i, sep = "")) 
      
    }
    
  } else if ( method == "DS" ) {
    for (i in 2:ncol(dinucleotide_table)) {
      
      colnames(dinucleotide_table)[i] <- c(paste(i - 8, "-", i - 7, sep = ""))
      
      x_order <- c(x_order, paste(i - 8, "-", i - 7, sep = "")) 
      
    }
  }
  dt_organized <- dinucleotide_table %>% gather(Position, count, 
                                                2:ncol(dinucleotide_table))
  
  # reorganize
  colnames(dt_organized) <- c("dinucleotides", "positions", "counts")
  
  dt_organized$freq = 100*dt_organized$counts/sum(dinucleotide_table[2])
  
  dt_organized <- filter(dt_organized, dinucleotides == "TT" |
                           dinucleotides == "TC" | dinucleotides == "CT" |
                           dinucleotides == "CC")
  
  dt_organized$dinucleotides <- factor(dt_organized$dinucleotides, 
                                       levels = c("CC", "CT", "TC", "TT"))
  
  dt_organized$positions <- factor(dt_organized$positions, 
                                   levels = x_order)
  
  return(dt_organized)
}


#### Main ####

flog.info("Importing length dist sample (plot B)...")
pB_data <- read.delim( argv$len_xr_cpd_120, header = FALSE )
colnames(pB_data) <- c("oligomer_length", "counts")

flog.info("Importing dinucleotide tables (plot C)...")
pC1_sample <-  read.table( argv$dinuc_xr_cpd_120, header = TRUE )
pC2_sample <- read.table( argv$dinuc_ds_cpd_120, header = TRUE )

flog.info("Reshaping dinucleotide tables (plot C)...")
pC1_data <- dinuc(pC1_sample, "XR")
pC2_data <- dinuc(pC2_sample, "DS") 

pC1_data$method <- "XR-seq"
pC1_data$product <- "CPD"
pC2_data$product <- "CPD"

flog.info("Importing file of read counts on TSS region (plot D)...")
pD_sample_df <- read.delim(argv$tss, header = FALSE, sep = "\t")

colnames(pD_sample_df) <- c("chromosomes", "start_position", "end_position", 
                            "dataset", "score", "dataset_strand", "counts", 
                            "sample_names", "file_names", "layout", "cell_line", 
                            "product", "method", "uv_exposure", "treatment", 
                            "phase", "time_after_exposure", "replicate", 
                            "project", "sample_source", "sample_strand", 
                            "mapped_reads", "RPKM")

pD_df_rr <- repair_rate( pD_sample_df )
pD_df_rr_org <- window_numbering( pD_df_rr, 4, 101 )
pD_df_rr_org$dataset <- "TSS"

pD_df_rr_org$dataset_strand <- factor(
  pD_df_rr_org$dataset_strand, levels = c("TS","NTS"))

flog.info("Importing file of read counts on TES region (plot D)...")
pD2_sample_df <- read.delim(argv$tes, header = FALSE, sep = "\t")

colnames(pD2_sample_df) <- c("chromosomes", "start_position", "end_position", 
                             "dataset", "score", "dataset_strand", "counts", 
                             "sample_names", "file_names", "layout", "cell_line", 
                             "product", "method", "uv_exposure", "treatment", 
                             "phase", "time_after_exposure", "replicate", 
                             "project", "sample_source", "sample_strand", 
                             "mapped_reads", "RPKM")

pD2_df_rr <- repair_rate( pD2_sample_df )
pD2_df_rr_org <- window_numbering( pD2_df_rr, 4, 101 )
pD2_df_rr_org$dataset <- "TES"

pD2_df_rr_org$dataset_strand <- factor(
  pD2_df_rr_org$dataset_strand, levels = c("TS","NTS"))

pD_comb <- rbind(pD_df_rr_org, pD2_df_rr_org)
pD_comb$dataset <- factor(pD_comb$dataset, levels = c("TSS", "TES"))

pD_comb$sample <- "CPD\n\n120 min."

#### Filtering Samples ####

flog.info("Filtering samples (plot D)...")
pD_data <- filter(pD_df_rr_org, phase != "async", replicate == "_",
                  time_after_exposure == "120")

pD2_data <- filter(pD2_df_rr_org, phase != "async", replicate == "_",
                   time_after_exposure == "120")

pD_comb_data <- filter(pD_comb, phase != "async", replicate == "_",
                       time_after_exposure == "120", phase == "early")

flog.info("Plotting...")
#### Plot A ####

flog.info("Plotting A...")
# plot A will be experimental setup drawing
# we are creating an empty text for that part
p.A <- wrap_elements(grid::textGrob(''))

flog.info("Plotting B...")
#### Plot B.1 ####

pB_data$sample <- "CPD 120 min."

write.table(pB_data, file = paste0(argv$data_prefix, "B.csv"), quote = FALSE, 
            row.names = FALSE, sep = ",")

pB_data$sample <- "CPD\n\n120 min."

# create the plot 
p.B <- ggplot(pB_data, aes(x = oligomer_length, y = counts/1000000)) + 
  #fill = counts)) + 
  geom_bar(stat = "identity") +
  facet_grid(sample~.) +
  scale_y_continuous(breaks = c(0, 6, 12),
                     labels = c("0", "6", expression(paste("(x10"^"6",")", 
                                                           sep = "")))) +
  xlim(15, 35) +
  xlab("Oligomer Length") + ylab("Counts") 

# adding and overriding the default plot format
p.B <- p.B + p_format 

flog.info("Plotting C...")
#### Plot C.1 ####

write.table(pC1_data, file = paste0(argv$data_prefix, "C1.csv"), 
            quote = FALSE, row.names = FALSE, sep = ",")

# create the plot 
p.C.1 <- ggplot(pC1_data, 
                aes(x = positions, y = freq, fill = dinucleotides)) + 
  geom_bar(stat = "identity") +
  facet_grid(~method) +
  scale_y_continuous(breaks = c(0, 50, 100),
                     limits = c(0, 101)) +
  scale_x_discrete(breaks = c("10-11","20-21"),
                   labels = c("10-11","20-21")) + 
  scale_fill_viridis(discrete = TRUE, begin = 1, end = 0) +
  ylab("Frequency (%)") + 
  labs(fill = "")

# adding and overriding the default plot format
p.C.1 <- p.C.1 + p_format +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(size = 12, margin=margin(0,0,-40,0)),
        legend.position = c(0.2, 0.9),
        legend.background = element_blank(),
        legend.key.size = unit(3, "mm")) 


#### Plot C.2 ####

pC2_data$method <- "Damage-seq"

write.table(pC2_data, file = paste0(argv$data_prefix, "C2.csv"), 
            quote = FALSE, row.names = FALSE, sep = ",")

pC2_data$method <- "Damage-\nseq"

# create the plot 
p.C.2 <- ggplot(pC2_data, 
                aes(x = positions, y = freq, fill = dinucleotides)) + 
  geom_bar(stat = "identity") +
  facet_grid(~product~method) +
  scale_y_continuous(breaks = c(0, 50, 100),
                     limits = c(0, 101)) +
  scale_x_discrete(breaks = c("-2--1"),
                   labels = c("-1")) + 
  scale_fill_viridis(discrete = TRUE, begin = 1, end = 0) +
  ylab("Frequency (%)") + 
  labs(fill = "")

# adding and overriding the default plot format
p.C.2 <- p.C.2 + p_format + guides(fill = "none") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 12, margin=margin(0,0,-40,0)),
        axis.title.x = element_blank()) 

flog.info("Plotting D...")
#### Plot D ####

# create the plot 
p.D <- ggplot(pD_data, aes(x = windows, y = log2(xr_ds))) + 
  geom_vline(xintercept = 0, color = "gray", linetype = "dashed") +
  geom_line(aes(color = dataset_strand, alpha = .9)) + 
  facet_grid(~product~time_after_exposure~phase,
             labeller = labeller(product = product_labs, 
                                 method = method_labs, 
                                 time_after_exposure = taex_labs, 
                                 replicate = rep_labs, phase = phase_labs)) + 
  xlab("Position Relative to TSS (kb)") + 
  ylab(fr_xr_ds_lab) +
  scale_y_continuous(breaks = c(0, 1, 2, 3),
                     limits = c(0, 3)) +
  scale_x_continuous(limits = c(-101, 101), 
                     breaks = c(-101, 0, 101), 
                     labels = c("-10", "0", "+10")) + 
  scale_color_manual(values = c("magenta2", "seagreen")) +
  labs(color = "") +
  guides(alpha = "none")

# create the plot 
p.D2 <- ggplot(pD2_data, aes(x = windows, y = log2(xr_ds))) + 
  geom_vline(xintercept = 0, color = "gray", linetype = "dashed") +
  geom_line(aes(color = dataset_strand, alpha = .9)) + 
  facet_grid(~product~time_after_exposure~phase,
             labeller = labeller(product = product_labs, 
                                 method = method_labs, 
                                 time_after_exposure = taex_labs, 
                                 replicate = rep_labs, phase = phase_labs)) + 
  xlab("Position Relative to TES (kb)") + 
  ylab(fr_xr_ds_lab) +
  scale_color_manual(values = c("magenta2", "seagreen")) +
  labs(color = "") +
  guides(alpha = "none")

pD_comb_data$sample <- "CPD 120 min."

write.table(pD_comb_data, file = paste0(argv$data_prefix, "D.csv"), 
            quote = FALSE, row.names = FALSE, sep = ",")

pD_comb_data$sample <- "CPD\n\n120 min."

# create the plot 
p.D_comb <- ggplot(pD_comb_data, aes(x = windows, y = log2(xr_ds))) + 
  geom_vline(xintercept = 0, color = "gray", linetype = "dashed") +
  geom_line(aes(color = dataset_strand, alpha = .9)) + 
  facet_grid(~sample~dataset) + 
  xlab("Position Relative to TSS and TES (kb)") + 
  scale_x_continuous(limits = c(-101, 101), 
                     breaks = c(-101, 0, 101), 
                     labels = c("-10", "0", "+10")) +  
  ylab(fr_xr_ds_lab) +
  scale_color_manual(values = c("magenta2", "seagreen")) +
  labs(color = "") +
  guides(alpha = "none")

# adding and overriding the default plot format
p.D_comb <- p.D_comb + p_format +
  theme(legend.position = c(0.5, 0.9),
        axis.text.x = element_text(hjust=c(0.1, 0.4, 0.5, 0.6, 0.9)))

flog.info("Plotting all together...")
#### Combining Plots with Patchwork ####

layout <- "
AAAABB
AAAABB
AAAABB
AAAABB
CCCDDD
CCCDDD
CCCDDD
CCCDDD
CCCDDD
"
layout2 <- "
AAAAAAAABBBB
"

p.A.B <- p.A + p.B

p.C <- p.C.1 + (p.C.2 + plot_layout(tag_level = 'new')) +
  plot_layout(design = layout2)

p_final <- p.A.B + p.C + p.D_comb +
  plot_layout(design = layout) +
  plot_annotation(tag_levels = 'A', caption = 
                    'Position Relative to the first base of reads',
                  theme = theme(plot.caption = element_text(size = 12,
                                                            hjust = .1, 
                                                            vjust = 9.5))) &
  theme(plot.tag = element_text(size = 12, face="bold"))

ggsave(argv$o, width = 22, height = 18, units = "cm")

