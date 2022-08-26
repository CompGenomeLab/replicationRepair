#### Packages and Libraries ####

library(stringr)
library(ggpubr)
library(ggplot2)
library(dplyr)
library(reshape2)
library(patchwork)
library(grid)
library(argparser)

######## Arguments ##########
p <- arg_parser("producing the figure 4C, 4D, S13, S15, S16B, and S16C")
p <- add_argument(p, "--df", help="region file with read counts")
p <- add_argument(p, "--df_sim", help="region file with simulated read counts")
p <- add_argument(p, "--phase", help="early or late")
p <- add_argument(p, "--intergenic", help="True if reads come from intergenic regions")
p <- add_argument(p, "--data_prefix", help="name prefix of the dataframes that generate the plots")
p <- add_argument(p, "-o", help="output")

# Parse the command line arguments
argv <- parse_args(p)

sample_csv <- argv$df

sample_sim_csv <- argv$df_sim

myphase <- argv$phase

if (argv$intergenic == "True"){ xlabname = "Position Relative to Initiation Zones (kb)\n(Intergenic)" 
} else if (argv$intergenic == "False"){ xlabname = "Position Relative to Initiation Zones (kb)"
}

if (myphase == "early"){ 
  ylabname1 = 'CPD\n12 min.\nEarly S Phase'
  ylabname2 = 'CPD\n120 min.\nEarly S Phase'
} else if (myphase == "late"){  
  ylabname1 = 'CPD\n12 min.\nLate S Phase'
  ylabname2 = 'CPD\n120 min.\nLate S Phase'
  }


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

df_rr <- repair_rate( sample_df )
df_rr_org <- window_numbering( df_rr, 4, 101 )
df_rr_org <- domain_name( df_rr_org, 1 )
df_rr_org$dataset <- "Initiation Zones"

df_rr_org$sample_strand <- factor(
  df_rr_org$sample_strand, levels = c("+","-"))

sample_df_sim <- read.delim( sample_sim_csv, header = F )

colnames(sample_df_sim) <- c("chromosomes", "start_position", "end_position", 
                         "dataset", "score", "dataset_strand", "counts", 
                         "sample_names", "file_names", "layout", "cell_line", 
                         "product", "method", "uv_exposure", "treatment", 
                         "phase", "time_after_exposure", "replicate", 
                         "project", "sample_source", "sample_strand", 
                         "mapped_reads", "RPKM")

df_rr_sim <- repair_rate( sample_df_sim )
df_rr_org_sim <- window_numbering( df_rr_sim, 4, 101 )
df_rr_org_sim <- domain_name( df_rr_org_sim, 1 )
df_rr_org_sim$dataset <- "Initiation Zones"

df_rr_org_sim$sample_strand <- factor(
  df_rr_org_sim$sample_strand, levels = c("+","-"))

colnames(df_rr_org)[12] <- "real"
colnames(df_rr_org_sim)[12] <- "sim"

df_rr_rs <- merge(df_rr_org, df_rr_org_sim, by = c("chromosomes", "start_position", 
                                           "end_position", "dataset", "score",
                                           "dataset_strand", "product", "windows", 
                                           "phase", "time_after_exposure", "repdomains", 
                                           "sample_strand", "replicate" ))
df_rr_rs$xr_ds <- df_rr_rs$real / df_rr_rs$sim

# filtering for B.1
pB1_data <- filter(df_rr_rs, phase != "async", replicate == "_", 
                   product == "CPD", time_after_exposure == "12", 
                   phase == myphase)


# filtering for C.1
pC1_data <- filter(df_rr_rs, phase != "async", replicate == "_", 
                   product == "CPD", time_after_exposure == "120", 
                   phase == myphase)


#### Plot A ####

# plot A will be a drawing
# we are creating an empty text for that part
p.A <- wrap_elements(grid::textGrob(''))


#### Plot B.1 ####

# create the plot 
p.B.1 <- ggplot(pB1_data, aes(x = windows, y = log2(xr_ds))) +
  geom_vline(xintercept = 0, color = "gray", linetype = "dashed") + 
  geom_line(aes(color = sample_strand)) + 
  facet_grid(~repdomains) +
  xlab(windows_lab) + ylab("n. Repair Rates (log2)") +
  scale_x_continuous(limits = c(-101, 101), 
                     breaks = c(-101, 0, 101), 
                     labels = c("-10", "0", "+10")) + 
  scale_color_manual(values = strand_colors) + 
  labs(color = "Strands")

# adding and overriding the default plot format
p.B.1 <- p.B.1 + p_format + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_text(hjust=c(0.3, 0.5, 0.7)))


#### Plot B.2 ####

pB1_data$direction[pB1_data$windows<0] <- "Left Replicating" 
pB1_data$direction[pB1_data$windows>0] <- "Right Replicating" 

pB1_data <- filter(pB1_data, direction != "NA")
pB1_data$log2val <- log2(pB1_data$xr_ds)

pB1_boxplot <- pB1_data[,c("repdomains", "log2val", "sample_strand", "direction", "windows")]

p.B.2 <- ggplot(data=pB1_boxplot, aes(x=repdomains, y=log2val, fill=sample_strand)) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~direction) +
  xlab("") + 
  ylab("") +
  scale_fill_manual(values = c("+" = "#0571b0", 
                               "-" = "#ca0020"), 
                    guide = "none") +
  #scale_y_continuous(breaks = c(0, 1, 2),
  #                   limits = c(0, 2)) +
  labs(color = "Strands", fill = "") 

# adding and overriding the default plot format
p.B.2 <- p.B.2 + p_format + 
  theme(strip.background = element_blank(),
        ) 


#### Plot C.1 ####

# create the plot 
p.C.1 <- ggplot(pC1_data, aes(x = windows, y = log2(xr_ds))) + 
  geom_vline(xintercept = 0, color = "gray", linetype = "dashed") +
  geom_line(aes(color = sample_strand)) + 
  facet_grid(~repdomains) +
  xlab(xlabname) + 
  ylab("n. Repair Rates (log2)") +
  scale_x_continuous(limits = c(-101, 101), 
                     breaks = c(-101, 0, 101), 
                     labels = c("-10", "0", "+10")) + 
  scale_color_manual(values = strand_colors) + 
  labs(color = "Strands")

# adding and overriding the default plot format
p.C.1 <- p.C.1 + p_format + 
  theme(axis.text.x=element_text(hjust=c(0.3, 0.5, 0.7)))


#### Plot C.2 ####

pC1_data$direction[pC1_data$windows<0] <- "Left Replicating" 
pC1_data$direction[pC1_data$windows>0] <- "Right Replicating" 

pC1_data <- filter(pC1_data, direction != "NA")
pC1_data$log2val <- log2(pC1_data$xr_ds)

pC1_boxplot <- pC1_data[,c("repdomains", "log2val", "sample_strand", "direction", "windows")]

p.C.2 <- ggplot(data=pC1_boxplot, aes(x=repdomains, y=log2val, fill=sample_strand)) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~direction) +
  xlab("Replication Domains") + 
  ylab("") +
  scale_fill_manual(values = c("+" = "#0571b0", 
                               "-" = "#ca0020"), 
                    guide = "none") +
  #scale_y_continuous(breaks = c(0, 1, 2),
  #                   limits = c(0, 2)) +
  labs(color = "Strands", fill = "") 

# adding and overriding the default plot format
p.C.2 <- p.C.2 + p_format + 
  theme(strip.background = element_blank(),
#        strip.text.x = element_blank(), 
#  theme(axis.title.x=element_blank(),
#        axis.text.x=element_blank(),
#        axis.ticks.x = element_blank(),
  )


#### Combining Plots with Patchwork ####

layout <- "
AAAAAAA
BBBCCCD
FFFGGGH
"
layout2 <- "
BBBCCCD
FFFGGGH
"

if (argv$intergenic == "False" & myphase == "late"){

write.table(pB1_data, file = paste0(argv$data_prefix, "C1.csv"), quote = FALSE, 
            row.names = FALSE, sep = ",")
write.table(pB1_boxplot, file = paste0(argv$data_prefix, "C2.csv"), quote = FALSE, 
            row.names = FALSE, sep = ",")
write.table(pC1_data, file = paste0(argv$data_prefix, "D1.csv"), quote = FALSE, 
            row.names = FALSE, sep = ",")
write.table(pC1_boxplot, file = paste0(argv$data_prefix, "D2.csv"), quote = FALSE, 
            row.names = FALSE, sep = ",")

p.B.2 <- p.B.2 + stat_compare_means(label = "p.signif",  paired = TRUE, label.y = 0.4) 
p.C.2 <- p.C.2 + stat_compare_means(label = "p.signif",  paired = TRUE, label.y = 0.05) 

(p.A + labs(title="A")) + 
(p.B.1 + labs(title="C")) + p.B.2 + 
  grid::textGrob(ylabname1, 
                  rot = -90, gp=gpar(fontsize=12), 
                  y = unit(.55, "npc")) + 
  (p.C.1 + labs(title="D")) + p.C.2 + 
  grid::textGrob(ylabname2, 
                  rot = -90, gp=gpar(fontsize=12), 
                  y = unit(.62, "npc")) + 
  plot_layout(design = layout, guides = "collect") & 
  theme(plot.tag = element_text(size = 12, face="bold"),
        legend.position = 'bottom', 
        plot.title = element_text(hjust = -0.2, 
                                  size = 12, face="bold"))

ggsave( argv$o, width = 22, height = 18, units = "cm" )

} else if (argv$intergenic == "True" & myphase == "early"){

write.table(pB1_data, file = paste0(argv$data_prefix, "B1.csv"), quote = FALSE, 
            row.names = FALSE, sep = ",")
write.table(pB1_boxplot, file = paste0(argv$data_prefix, "B2.csv"), quote = FALSE, 
            row.names = FALSE, sep = ",")
write.table(pC1_data, file = paste0(argv$data_prefix, "C1.csv"), quote = FALSE, 
            row.names = FALSE, sep = ",")
write.table(pC1_boxplot, file = paste0(argv$data_prefix, "C2.csv"), quote = FALSE, 
            row.names = FALSE, sep = ",")

p.B.2 <- p.B.2 + stat_compare_means(label = "p.signif",  paired = TRUE, label.y = 0.85) 
p.C.2 <- p.C.2 + stat_compare_means(label = "p.signif",  paired = TRUE, label.y = 0.15) 

p.A + labs(title="A") + 
(p.B.1 + labs(title="B")) + p.B.2 + 
  grid::textGrob(ylabname1, 
                  rot = -90, gp=gpar(fontsize=12), 
                  y = unit(.55, "npc")) + 
  (p.C.1 + labs(title="C")) + p.C.2 + 
  grid::textGrob(ylabname2, 
                  rot = -90, gp=gpar(fontsize=12), 
                  y = unit(.62, "npc")) + 
  plot_layout(design = layout, guides = "collect") & 
  theme(plot.tag = element_text(size = 12, face="bold"),
        legend.position = 'bottom', 
        plot.title = element_text(hjust = -0.2, 
                                  size = 12, face="bold"))

ggsave( argv$o, width = 22, height = 18, units = "cm" )

} else {

write.table(pB1_data, file = paste0(argv$data_prefix, "A1.csv"), quote = FALSE, 
            row.names = FALSE, sep = ",")
write.table(pB1_boxplot, file = paste0(argv$data_prefix, "A2.csv"), quote = FALSE, 
            row.names = FALSE, sep = ",")
write.table(pC1_data, file = paste0(argv$data_prefix, "B1.csv"), quote = FALSE, 
            row.names = FALSE, sep = ",")
write.table(pC1_boxplot, file = paste0(argv$data_prefix, "B2.csv"), quote = FALSE, 
            row.names = FALSE, sep = ",")

p.B.2 <- p.B.2 + stat_compare_means(label = "p.signif",  paired = TRUE) 
p.C.2 <- p.C.2 + stat_compare_means(label = "p.signif",  paired = TRUE) 

(p.B.1 + labs(title="A")) + p.B.2 + 
  grid::textGrob(ylabname1, 
                  rot = -90, gp=gpar(fontsize=12), 
                  y = unit(.55, "npc")) + 
  (p.C.1 + labs(title="B")) + p.C.2 + 
  grid::textGrob(ylabname2, 
                  rot = -90, gp=gpar(fontsize=12), 
                  y = unit(.62, "npc")) + 
  plot_layout(design = layout2, guides = "collect") & 
  theme(plot.tag = element_text(size = 12, face="bold"),
        legend.position = 'bottom', 
        plot.title = element_text(hjust = -0.2, 
                                  size = 12, face="bold"))

ggsave( argv$o, width = 22, height = 18, units = "cm" )

}