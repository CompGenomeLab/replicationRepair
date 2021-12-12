#### Packages and Libraries ####

library(argparser)
library(stringr)
library(ggpubr)
library(ggplot2)
library(dplyr)
library(reshape2)
library(patchwork)
library(grid)

######## Arguments ##########
p <- arg_parser("producing the suplementary figure 5")
p <- add_argument(p, "--real", help="windowed (20kb) initiation zones file with read counts")
p <- add_argument(p, "--sim", help="windowed (20kb) initiation zones file with simulated read counts")
p <- add_argument(p, "--prod", help="Damage type (64_PP or CPD)")
p <- add_argument(p, "-o", help="output")


# Parse the command line arguments
argv <- parse_args(p)

#### Variables ####

rep <- "_" 

prod <- argv$prod

#### Default Plot Format ####

source("workflow/scripts/plot_format.R")


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
pA1_data <- filter(real_df_org, method != "DNA_seq", phase != "async", 
                   replicate == rep, method == "Damage_seq", product == prod)

# filtering for A.2
pA2_data <- filter(sim_df_org, phase != "async", replicate == rep, 
                   product == prod, method == "Damage_seq")

# filtering for B.1
pB1_data <- filter(real_df_org, method != "DNA_seq", phase != "async", 
                   replicate == rep, method == "XR_seq", product == prod)

# filtering for B.2
pB2_data <- filter(sim_df_org, phase != "async", method != "DNA_seq", 
                   replicate == rep, product == prod, method == "XR_seq")

# for naming of simulated samples
method_labs_sim <- c("Simulated \nDamage-seq", "Simulated \nXR-seq")
names(method_labs_sim) <- c("Damage_seq", "XR_seq")

#### Plot A.1 ####

# create the plot 
p.A.1 <- ggplot(pA1_data, aes(x = windows, y = RPKM)) + 
  geom_vline(xintercept = 0, color = "gray", linetype = "dashed") +
  geom_line(aes(color = sample_strand)) + 
  facet_grid(~repdomains~product~time_after_exposure~phase~method, 
             labeller = labeller(product = product_labs, 
                                 method = method_labs, 
                                 time_after_exposure = taex_labs,
                                 phase = phase_labs)) + 
  xlab("") + ylab(fr_lab) +
  scale_y_continuous(breaks = c(.1, .3, .5),
                     limits = c(.0, .5)) +
  scale_x_continuous(limits = c(-101, 101), 
                     breaks = c(-101, 0, 101), 
                     labels = c("-10", "0", "+10")) + 
  scale_color_manual(values = strand_colors) + 
  labs(color = "Strands")

# adding and overriding the default plot format
p.A.1 <- p.A.1 + p_format +
  theme(strip.text.y = element_text(size=0, margin = margin(0, 0, 0, 0)),
        panel.border = element_rect(fill = NA)) 


#### Plot A.2 ####

# create the plot
p.A.2 <- ggplot(pA2_data, aes(x = windows, y = RPKM)) + 
  geom_line(aes(color = sample_strand)) + 
  geom_vline(xintercept = 0, color = "gray", linetype = "dashed") +
  facet_grid(~repdomains~product~time_after_exposure~phase~method,
             labeller = labeller(product = product_labs,
                                 method = method_labs_sim,
                                 time_after_exposure = taex_labs,
                                 phase = phase_labs)) +
  xlab("") + ylab(fr_lab) +
  scale_y_continuous(breaks = c(.1, .3, .5),
                     limits = c(.0, .5)) +
  scale_x_continuous(limits = c(-101, 101), 
                     breaks = c(-101, 0, 101), 
                     labels = c("-10", "0", "+10")) + 
  scale_color_manual(values = strand_colors) + 
  labs(color = "Strands")

# adding and overriding the default plot format
p.A.2 <- p.A.2 + p_format + 
  theme(axis.title.y=element_blank(),
        panel.border = element_rect(fill = NA),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        strip.text.y = element_text(size=0, margin = margin(0, 0, 0, 0)))


#### Plot B.1 ####

# create the plot 
p.B.1 <- ggplot(pB1_data, aes(x = windows, y = RPKM)) + 
  geom_vline(xintercept = 0, color = "gray", linetype = "dashed") +
  geom_line(aes(color = sample_strand)) + 
  facet_grid(~repdomains~product~time_after_exposure~phase~method, 
             labeller = labeller(product = product_labs, 
                                 method = method_labs, 
                                 time_after_exposure = taex_labs,
                                 phase = phase_labs)) + 
  xlab("") + ylab(fr_lab) +
  scale_y_continuous(breaks = c(.1, .5, 1),
                     limits = c(.0, 1)) +
  scale_x_continuous(limits = c(-101, 101), 
                     breaks = c(-101, 0, 101), 
                     labels = c("-10", "0", "+10")) + 
  scale_color_manual(values = strand_colors) + 
  labs(color = "Strands")

# adding and overriding the default plot format
p.B.1 <- p.B.1 + p_format + 
  theme(strip.text.y = element_text(size=0, margin = margin(0, 0, 0, 0)),
        axis.title.y = element_blank(),
        panel.border = element_rect(fill = NA)) 


#### Plot B.2 ####

# create the plot 
p.B.2 <- ggplot(pB2_data, aes(x = windows, y = RPKM)) + 
  geom_line(aes(color = sample_strand)) + 
  geom_vline(xintercept = 0, color = "gray", linetype = "dashed") +
  facet_grid(~repdomains~product~time_after_exposure~phase~method, 
             labeller = labeller(product = product_labs, 
                                 method = method_labs_sim, 
                                 time_after_exposure = taex_labs, 
                                 phase = phase_labs)) + 
  xlab("") + ylab(fr_lab) +
  #scale_y_continuous(breaks = c(.1, .5, 1),
  #                   limits = c(.0, 1)) +
  scale_x_continuous(limits = c(-101, 101), 
                     breaks = c(-101, 0, 101), 
                     labels = c("-10", "0", "+10")) + 
  scale_color_manual(values = strand_colors) + 
  labs(color = "Strands")

# adding and overriding the default plot format
#p.B.2 <- p.B.2 + p_format + 
#  theme(axis.title.y=element_blank(),
#        axis.text.y=element_blank(),
#        axis.ticks.y=element_blank(),
#        panel.border = element_rect(fill = NA))


#### Combining Plots with Patchwork ####

p.A.2 <- p.A.2 + plot_layout(tag_level = 'new') 
p.B.2 <- p.B.2 + plot_layout(tag_level = 'new') 

(p.A.1 | p.A.2 | p.B.1 | p.B.2) + 
  plot_annotation(caption = 
                    'Position Relative to Initiation Zones (kb)',
                  theme = theme(plot.caption = 
                                  element_text(size = 12, 
                                               hjust = .43, vjust = 20))) +
  plot_layout(guides = "collect") & 
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size = 12, face="bold"),
        legend.position = 'bottom', 
        plot.title = element_text(hjust = -0.2, vjust = 5, 
                                  size = 12, face="bold"))


ggsave(argv$o, width = 22, height = 18, units = "cm")

