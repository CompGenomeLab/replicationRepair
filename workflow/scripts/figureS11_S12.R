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
p <- arg_parser("producing the figure S8 and S9")
p <- add_argument(p, "--df", help="region file with read counts")
p <- add_argument(p, "--df_sim", help="region file with simulated read counts")
p <- add_argument(p, "--intergenic", help="True if reads come from intergenic regions")
p <- add_argument(p, "--data_prefix", help="name prefix of the dataframes that generate the plots")
p <- add_argument(p, "-o", help="output")

# Parse the command line arguments
argv <- parse_args(p)

sample_csv <- argv$df

sample_sim_csv <- argv$df_sim

if (argv$intergenic == "True"){ xlabname = "Position Relative to Initiation Zones (kb)\n(Intergenic)" 
} else if (argv$intergenic == "False"){ xlabname = "Position Relative to Initiation Zones (kb)"
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

# filtering for A.1
pA1_data <- filter(df_rr_rs, phase != "async", replicate == "_", 
                   product == "64_PP", phase == "early")

# for plot A.2
pA2_data <- rr_boxplot( pA1_data ) 


# filtering for B.1
pB1_data <- filter(df_rr_rs, phase != "async", replicate == "_", 
                   product == "64_PP", phase == "late")

# for plot B.2
pB2_data <- rr_boxplot( pB1_data ) 


#### Plot A.1 ####

write.table(pA1_data, file = paste0(argv$data_prefix, "A1.csv"), 
            quote = FALSE, row.names = FALSE, sep = ",")

# create the plot 
p.A.1 <- ggplot(pA1_data, aes(x = windows, y = log2(xr_ds))) +
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
p.A.1 <- p.A.1 + p_format + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_text(hjust=c(0.3, 0.5, 0.7)))


#### Plot A.2 ####

pA1_data$direction[pA1_data$windows<0] <- "Left Replicating" 
pA1_data$direction[pA1_data$windows>0] <- "Right Replicating" 

pA1_data <- filter(pA1_data, direction != "NA")
pA1_data$log2val <- log2(pA1_data$xr_ds)

pA1_boxplot <- pA1_data[,c("repdomains", "log2val", "sample_strand", "direction", "windows")]

write.table(pA1_boxplot, file = paste0(argv$data_prefix, "A2.csv"), 
            quote = FALSE, row.names = FALSE, sep = ",")

p.A.2 <- ggplot(data=pA1_boxplot, aes(x=repdomains, y=log2val, fill=sample_strand)) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~direction) +
  xlab("") + 
  ylab("") +
  scale_fill_manual(values = c("+" = "#0571b0", 
                               "-" = "#ca0020"), 
                    guide = "none") +
  labs(color = "Strands", fill = "") 

# adding and overriding the default plot format
p.A.2 <- p.A.2 + p_format + 
  stat_compare_means(label = "p.signif",  paired = TRUE) + 
  theme(strip.background = element_blank()) 


#### Plot B.1 ####

write.table(pB1_data, file = paste0(argv$data_prefix, "B1.csv"), 
            quote = FALSE, row.names = FALSE, sep = ",")

# create the plot 
p.B.1 <- ggplot(pB1_data, aes(x = windows, y = log2(xr_ds))) + 
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
p.B.1 <- p.B.1 + p_format + 
  theme(axis.text.x=element_text(hjust=c(0.3, 0.5, 0.7)))


#### Plot B.2 ####

pB1_data$direction[pB1_data$windows<0] <- "Left Replicating" 
pB1_data$direction[pB1_data$windows>0] <- "Right Replicating" 

pB1_data <- filter(pB1_data, direction != "NA")
pB1_data$log2val <- log2(pB1_data$xr_ds)

pB1_boxplot <- pB1_data[,c("repdomains", "log2val", "sample_strand", "direction", "windows")]

write.table(pB1_boxplot, file = paste0(argv$data_prefix, "B2.csv"), 
            quote = FALSE, row.names = FALSE, sep = ",")

p.B.2 <- ggplot(data=pB1_boxplot, aes(x=repdomains, y=log2val, fill=sample_strand)) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~direction) +
  xlab("Replication Domains") + 
  ylab("") +
  scale_fill_manual(values = c("+" = "#0571b0", 
                               "-" = "#ca0020"), 
                    guide = "none") +
  labs(color = "Strands", fill = "") 

# adding and overriding the default plot format
p.B.2 <- p.B.2 + p_format + 
  stat_compare_means(label = "p.signif",  paired = TRUE) + 
  theme(strip.background = element_blank())


#### Combining Plots with Patchwork ####

layout <- "
BBBCCCD
FFFGGGH
"
(p.A.1 + labs(title="A")) + p.A.2 + 
  grid::textGrob('(6-4)PP\n12 min.\nEarly S Phase', 
                  rot = -90, gp=gpar(fontsize=12), 
                  y = unit(.55, "npc")) + 
  (p.B.1 + labs(title="B")) + p.B.2 + 
    grid::textGrob('(6-4)PP\n12 min.\nLate S Phase', 
                    rot = -90, gp=gpar(fontsize=12), 
                    y = unit(.62, "npc")) + 
  plot_layout(design = layout, guides = "collect") & 
  theme(legend.position = 'bottom', 
        plot.title = element_text(hjust = -0.2, 
                                  size = 12, face="bold"))

ggsave( argv$o, width = 22, height = 18, units = "cm" )
