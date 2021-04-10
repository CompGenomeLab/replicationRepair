#### Packages and Libraries ####

library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(grid)
library(gridExtra)
library("viridis")
library(scales)


#### Variables ####

# name of the sample bed files for plot B
dinuc_ds_cpd_120 <- paste("~/Documents/myprojects/replicationRepair/3_output/",
                        "gitignore/2_sample_control/dinucleotide_composition/", 
                        "XT4H1-CT5H2-HelaD3-5R2h2_combined_cutadapt_sorted_",
                        "10_dinucleotideTable.txt", sep = "")

dinuc_xr_cpd_120 <- paste("~/Documents/myprojects/replicationRepair/3_output/",
                         "gitignore/2_sample_control/dinucleotide_composition/", 
                         "R19026421-2019XR1-Hela35X3_combined_R1_cutadapt_",
                         "sorted_26_dinucleotideTable.txt", sep = "")

# name of the sample txt files for plot C
len_xr_cpd_120 <- paste("~/Documents/myprojects/replicationRepair/3_output/", 
                       "gitignore/2_sample_control/length_distribution/", 
                       "R19026421-2019XR1-Hela35X3_combined_R1_cutadapt_",
                       "length_distribution.txt", 
                       sep = "")

# name of the sample csv file for plot D  
sample_csv_pD <- paste("~/Documents/myprojects/replicationRepair/3_output/",
                       "gitignore/1_TextforPlotting/", 
                       "[2020.11.25]final_report_tss_ready.csv", 
                       sep = "")

# path of the default plot format and functions
sourcePath <- "~/Documents/myprojects/replicationRepair/1_code/r/"


#### Default Plot Format ####

source("4_plot_format.R")


#### Fuctions ####

source("4_functions.R")

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

# for plot B
pB_data <- read.delim( len_xr_cpd_120, header = FALSE )
colnames(pB_data) <- c("oligomer_length", "counts")
pB_data$sample <- "CPD\n\n120 min."

# for plot C
pC1_sample <-  read.table( dinuc_xr_cpd_120, header = TRUE )
pC2_sample <- read.table( dinuc_ds_cpd_120, header = TRUE )

pC1_data <- dinuc(pC1_sample, "XR")
pC2_data <- dinuc(pC2_sample, "DS") 

pC1_data$method <- "XR-seq"
pC2_data$method <- "Damage-\nseq"
pC2_data$product <- "CPD"

# for plot D
pD_sample_df <- read.csv( sample_csv_pD )
pD_df_rr <- repair_rate( pD_sample_df )
pD_df_rr_org <- window_numbering( pD_df_rr, 4, 101 )
pD_df_rr_org$dataset <- "TSS"

pD_df_rr_org$dataset_strand <- factor(
  pD_df_rr_org$dataset_strand, levels = c("TS","NTS"))


#### Filtering Samples ####

pD_data <- filter(pD_df_rr_org, phase != "async", replicate == "A",
                  time_after_exposure == "120")

pD_data$sample <- "CPD\n\n120 min."


#### Plot A ####

# plot A will be experimental setup drawing
# we are creating an empty text for that part
p.A <- wrap_elements(grid::textGrob(''))


#### Plot B ####

# create the plot 
p.B <- ggplot(pB_data, aes(x = oligomer_length, y = counts/1000000)) +
                           #fill = counts)) + 
  geom_bar(stat = "identity") +
  facet_grid(sample~.) +
  scale_y_continuous(breaks = c(0, 5, 10, 15),
                     labels = c("0", "5", "10", expression(paste("(x10"^"6",")", 
                                                                 sep = ""))),
                     limits = c(0, 15)) +
  #scale_fill_gradient(name = "",
  #                    limits = c(0, 13000000), 
  #                    breaks = c(0, 5000000, 10000000, 13000000),
  #                    labels = c(0, 5, 10, expression(paste("(x10"^"6",")", 
  #                                                            sep = "")))) +
  xlim(15, 35) +
  xlab("Oligomer Length") + ylab("Counts") +
  labs(fill = "")

# adding and overriding the default plot format
p.B <- p.B + p_format +
  theme(legend.text = element_text(size = 8, vjust = 0.1, hjust = 0.1))


#### Plot C.1 ####

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
        legend.position = c(0.2, 0.7),
        plot.title = element_text(size = 12, margin=margin(0,0,-40,0)))


#### Plot C.2 ####

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
p.C.2 <- p.C.2 + p_format + guides(fill = FALSE) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.y = element_blank())


#### Plot D ####

# create the plot 
p.D <- ggplot(pD_data, aes(x = windows, y = log2(xr_ds))) + 
  geom_vline(xintercept = 0, color = "gray", linetype = "dashed") +
  geom_line(aes(color = dataset_strand, alpha = .9)) + 
  facet_grid(~sample~phase, labeller = labeller(phase = phase_labs)) + 
  xlab("Position Relative to TSS (kb)") + 
  ylab(fr_xr_ds_lab) +
  scale_y_continuous(breaks = c(0, 1, 2),
                     limits = c(0, 3)) +
  scale_x_continuous(limits = c(-101, 101), 
                     breaks = c(-101, 0, 101), 
                     labels = c("-10", "0", "+10")) + 
  scale_color_manual(values = c("magenta2", "seagreen")) +
  labs(color = "") +
  guides(alpha = FALSE)

# adding and overriding the default plot format
p.D <- p.D + p_format +
  theme(legend.position = c(0.5, 0.9),
        axis.text.x = element_text(hjust=c(0.1, 0.5, 0.9)))
p.D

#### Combining Plots with Patchwork ####

layout <- "
AAAAABB
AAAAABB
AAAAABB
AAAAABB
CCCCCCC
CCCCCCC
CCCCCCC
CCCCCCC
CCCCCCC
"

layout2 <- "
AAABCCCC
"

p.C.D <- p.C.1 + p.C.2 + p.D +
  plot_layout(design = layout2)


p.A + p.B + p.C.D +
  plot_layout(design = layout) & 
  plot_annotation(caption = 
                    'Position Relative to the first base of reads',
                  theme = theme(plot.caption = element_text(size = 12,
                                                            hjust = .15, 
                                                            vjust = 9.5))) &
  theme(plot.tag = element_text(size = 12, face="bold"))

ggsave("~/Desktop/supfig2.pdf", width = 22, height = 18, units = "cm")

