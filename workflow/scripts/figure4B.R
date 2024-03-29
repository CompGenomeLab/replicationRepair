#### Packages and Libraries ####

library(argparser)
library(stringr)
library(ggpubr)
library(ggplot2)
library(dplyr)
library(reshape2)
library(patchwork)
library(grid)
set.seed(1) 

######## Arguments ##########
p <- arg_parser("producing figure 4B")
p <- add_argument(p, "--real", help="windowed (20kb) initiation zones file with read counts")
p <- add_argument(p, "--sim", help="windowed (20kb) initiation zones file with simulated read counts")
p <- add_argument(p, "--prod", help="Damage type (64_PP or CPD)")
p <- add_argument(p, "--data_prefix", help="name prefix of the dataframes that generate the plots")
p <- add_argument(p, "-o", help="output")


# Parse the command line arguments
argv <- parse_args(p)

#### Variables ####

rep <- "_" 

prod <- argv$prod

#### Default Plot Format ####

source("workflow/scripts/plot_format.R")

taex_labs <- c("0 min.", "0/12\n min.", "120\n min.", "60 min.")
names(taex_labs) <- c("0", "12", "120", "60")

#### Functions ####

source("workflow/scripts/functions.R")


#### Main ####

# for plot A.1 and B.1
real_df <- read.delim( argv$real, header = F )
colnames(real_df) <- c("chromosomes", "start_position", "end_position", 
                       "dataset", "score", "dataset_strand", "counts", 
                       "sample_names", "file_names", "layout", "cell_line", 
                       "product", "method", "uv_exposure", "treatment", 
                       "phase", "time_after_exposure", "replicate", 
                       "project", "sample_source", "sample_strand", 
                       "mapped_reads", "RPKM")

real_df<- filter( real_df, method != "DNA_seq", phase != "async", 
                 replicate == rep, product == prod)

real_df_org <- window_numbering( real_df, 4, 101 )

real_df_org$repdomains <- data.frame(str_split_fixed(real_df_org$dataset, "_", -1))[,2]

real_df_org$repdomains <- factor(real_df_org$repdomains, levels = c("UTZ", "ERD", "DTZ", "LRD"))

real_df_org$dataset <- gsub("_.*", "", real_df_org$dataset)

real_df_org$sample_strand <- factor(
  real_df_org$sample_strand, levels = c("+","-"))

# for plot A.2 and B.2
sim_df <- read.delim( argv$sim, header = F )
colnames(sim_df) <- c("chromosomes", "start_position", "end_position", 
                      "dataset", "score", "dataset_strand", "counts", 
                      "sample_names", "file_names", "layout", "cell_line", 
                      "product", "method", "uv_exposure", "treatment", 
                      "phase", "time_after_exposure", "replicate", 
                      "project", "sample_source", "sample_strand", 
                      "mapped_reads", "RPKM")

sim_df<- filter( sim_df, method != "DNA_seq", phase != "async", 
                replicate == rep, product == prod)

sim_df_org <- window_numbering( sim_df, 4, 101 )

sim_df_org$repdomains <- data.frame(str_split_fixed(sim_df_org$dataset, "_", -1))[,2]

sim_df_org$repdomains <- factor(sim_df_org$repdomains, levels = c("UTZ", "ERD", "DTZ", "LRD"))

sim_df_org$dataset <- gsub("_.*", "", sim_df$dataset)

sim_df_org$sample_strand <- factor(
  sim_df_org$sample_strand, levels = c("+","-"))


# filtering for A.1
pB1_data <- filter(real_df_org, method != "DNA_seq", phase == "late", repdomains == "LRD", time_after_exposure == "12",
                   replicate == rep, method == "Damage_seq", product == prod)

# filtering for A.2
pB2_data <- filter(sim_df_org, phase == "late", replicate == rep, repdomains == "LRD", time_after_exposure == "12",
                   product == prod, method == "Damage_seq")

# filtering for B.1
pB3_data <- filter(real_df_org, method != "DNA_seq", phase == "late", repdomains == "LRD", time_after_exposure == "12",
                   replicate == rep, method == "XR_seq", product == prod)

# filtering for B.2
pB4_data <- filter(sim_df_org, phase == "late", method != "DNA_seq", repdomains == "LRD", time_after_exposure == "12",
                   replicate == rep, product == prod, method == "XR_seq")

# for naming of simulated samples
method_labs_sim <- c("Simulated \nDamage-seq", "Simulated \nXR-seq")
names(method_labs_sim) <- c("Damage_seq", "XR_seq")
method_labs <- c("Damage-\n seq", "XR-seq")
names(method_labs) <- c("Damage_seq", "XR_seq")

#### Plot B.1 ####

write.table(pB1_data, file = paste0(argv$data_prefix, "B1.csv"), quote = FALSE, 
            row.names = FALSE, sep = ",")

# create the plot 
p.B.1 <- ggplot(pB1_data, aes(x = windows, y = RPKM)) + 
  geom_vline(xintercept = 0, color = "gray", linetype = "dashed") +
  geom_line(aes(color = sample_strand)) + 
  facet_grid(~method, 
             labeller = labeller(product = product_labs, 
                                 method = method_labs, 
                                 time_after_exposure = taex_labs,
                                 phase = phase_labs)) + 
  xlab("") + ylab(fr_lab) +
  scale_y_continuous(breaks = c(.1, .2, .3),
                     limits = c(.0, .3)) +
  scale_x_continuous(limits = c(-101, 101), 
                     breaks = c(-101, 0, 101), 
                     labels = c("-10", "0", "+10")) + 
  scale_color_manual(values = strand_colors) + 
  labs(color = "Strands")

# adding and overriding the default plot format
p.B.1 <- p.B.1 + p_format + ggtitle("Real") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.margin = margin(0, 0, 0, 0, "pt"),
        strip.text.y = element_text(size=0, margin = margin(0, 0, 0, 0)),
        panel.border = element_rect(fill = NA),
        axis.text.x = element_text(size = 10, vjust = 0.6, hjust = c(0.1, 0.5, 0.9)),
        strip.text.x = element_text(size=0, margin = margin(0, 0, 0, 0))) 

#### Plot B.2 ####

write.table(pB2_data, file = paste0(argv$data_prefix, "B2.csv"), quote = FALSE, 
            row.names = FALSE, sep = ",")

# create the plot
p.B.2 <- ggplot(pB2_data, aes(x = windows, y = RPKM)) + 
  geom_line(aes(color = sample_strand)) + 
  geom_vline(xintercept = 0, color = "gray", linetype = "dashed") +
  facet_grid(method~., 
             labeller = labeller(product = product_labs, 
                                 method = method_labs, 
                                 time_after_exposure = taex_labs,
                                 phase = phase_labs)) + 
  xlab("") + ylab(fr_lab) +
  scale_y_continuous(breaks = c(.1, .2, .3),
                     limits = c(.0, .3)) +
  scale_x_continuous(limits = c(-101, 101), 
                     breaks = c(-101, 0, 101), 
                     labels = c("-10", "0", "+10")) + 
  scale_color_manual(values = strand_colors) + 
  labs(color = "Strands")

# adding and overriding the default plot format
p.B.2 <- p.B.2 + p_format + ggtitle("Simulated") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.margin = margin(0, 0, 0, 0, "pt"),
        axis.title.y=element_blank(),
        panel.border = element_rect(fill = NA),
        axis.text.y=element_blank(),
        axis.text.x = element_text(size = 10, vjust = 0.6, hjust = c(0.1, 0.5, 0.9)),
        axis.ticks.y=element_blank())


#### Plot B.3 ####

write.table(pB3_data, file = paste0(argv$data_prefix, "B3.csv"), quote = FALSE, 
            row.names = FALSE, sep = ",")

# create the plot 
p.B.3 <- ggplot(pB3_data, aes(x = windows, y = RPKM)) + 
  geom_vline(xintercept = 0, color = "gray", linetype = "dashed") +
  geom_line(aes(color = sample_strand)) + 
  facet_grid(~method, 
             labeller = labeller(product = product_labs, 
                                 method = method_labs, 
                                 time_after_exposure = taex_labs,
                                 phase = phase_labs)) + 
  xlab("") + ylab(fr_lab) +
  scale_y_continuous(breaks = c(.1, .2, .3),
                     limits = c(.0, .3)) +
  scale_x_continuous(limits = c(-101, 101), 
                     breaks = c(-101, 0, 101), 
                     labels = c("-10", "0", "+10")) + 
  scale_color_manual(values = strand_colors) + 
  labs(color = "Strands")

# adding and overriding the default plot format
p.B.3 <- p.B.3 + p_format + 
  theme(strip.text.y = element_text(size=0, margin = margin(0, 0, 0, 0)),
        panel.border = element_rect(fill = NA),
        plot.margin = margin(0, 0, 0, 0, "pt"),
        axis.text.x = element_text(size = 10, vjust = 0.6, hjust = c(0.1, 0.5, 0.9)),
        strip.text.x = element_text(size=0, margin = margin(0, 0, 0, 0))) 


#### Plot B.4 ####

write.table(pB4_data, file = paste0(argv$data_prefix, "B4.csv"), quote = FALSE, 
            row.names = FALSE, sep = ",")

# create the plot 
p.B.4 <- ggplot(pB4_data, aes(x = windows, y = RPKM)) + 
  geom_line(aes(color = sample_strand)) + 
  geom_vline(xintercept = 0, color = "gray", linetype = "dashed") +
  facet_grid(method~., 
             labeller = labeller(product = product_labs, 
                                 method = method_labs, 
                                 time_after_exposure = taex_labs,
                                 phase = phase_labs)) + 
  xlab("") + ylab(fr_lab) +
  scale_y_continuous(breaks = c(.1, .2, .3),
                     limits = c(.0, .3)) +
  scale_x_continuous(limits = c(-101, 101), 
                     breaks = c(-101, 0, 101), 
                     labels = c("-10", "0", "+10")) + 
  scale_color_manual(values = strand_colors) + 
  labs(color = "Strands")

# adding and overriding the default plot format
p.B.4 <- p.B.4 + p_format + 
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x = element_text(size = 10, vjust = 0.6, hjust = c(0.1, 0.5, 0.9)),
        axis.ticks.y=element_blank(), plot.margin = margin(0, 0, 0, 0, "pt"),
        panel.border = element_rect(fill = NA))

#### Combining Plots with Patchwork ####

(p.B.1 + p.B.2 ) / (p.B.3 + p.B.4) +
  plot_annotation(caption = 
                    'Position Relative to Initiation Zones (kb)',
                  theme = theme(plot.caption = 
                                  element_text(size = 12, 
                                               hjust = .43, vjust = 20))) +
  plot_layout(guides = "collect") & 
  theme(plot.tag = element_text(size = 12, face="bold"),
        legend.position = 'bottom')


ggsave(argv$o, width = 11, height = 9, units = "cm")