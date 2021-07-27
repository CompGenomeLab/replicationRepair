#### Packages and Libraries ####

library(stringr)
library(ggpubr)
library(ggplot2)
library(dplyr)
library(reshape2)
library(patchwork)
library(grid)
library(ggthemes)


#### Variables ####

rep <- "A" 
if (rep == "A"){ fig_name = "~/Desktop/supfig11.svg" 
} else if (rep == "B"){ fig_name = "no_need.svg" 
} 

# name of the sample csv file for plot A
real_csv <- paste("/home/azgarian/Documents/myprojects/replicationRepair/",
                  "final/final_reports_hg19_",
                  "iz_hela_windows_201_100_ready.csv", 
                  sep = "")

# name of the sample csv file for plot B
sim_csv_pA2 <- paste("~/Documents/myprojects/replicationRepair/3_output/",
                     "gitignore/1_TextforPlotting/", 
                     "[2020.01.19]final_report_InZones_windows_sim_ready.csv", 
                     sep = "")

sim_csv_pB2 <- paste("~/Documents/myprojects/replicationRepair/3_output/",
                     "gitignore/1_TextforPlotting/", 
                     "[2020.02.25]final_report_inZones_windows_201_100_",
                     "sim_ready.csv", sep = "")

sim_csv <- paste("/home/azgarian/Documents/myprojects/replicationRepair/",
                 "final/final_reports_sim_hg19_",
                 "iz_hela_windows_201_100_ready.csv", 
                 sep = "")

#### Default Plot Format ####

source("/home/azgarian/Documents/myprojects/replicationRepair/1_code/4_plot_format.R")


#### Fuctions ####

source("/home/azgarian/Documents/myprojects/replicationRepair/1_code/4_functions.R")


#### Main ####

# for plot A.1 and B.1
real_df <- read.csv( real_csv )

real_df<- filter(real_df, method != "DNA_seq", phase != "async", 
                 replicate == "A", product == "CPD")

real_df_org <- window_numbering( real_df, 4, 101 )
real_df_org$dataset <- gsub("_.*", "", real_df_org$dataset)

real_df_org$sample_strand <- factor(
  real_df_org$sample_strand, levels = c("+","-"))

# for plot A.2 and B.2
sim_df <- read.csv( sim_csv )

sim_df<- filter(sim_df, method != "DNA_seq", phase != "async", 
                 replicate == "A", product == "CPD")

sim_df_org <- window_numbering( sim_df, 4, 101 )
sim_df_org$dataset <- gsub("_.*", "", sim_df$dataset)

sim_df_org$sample_strand <- factor(
  sim_df_org$sample_strand, levels = c("+","-"))


# filtering for A.1
pA1_data <- filter(real_df_org, method != "DNA_seq", phase != "async", 
                   replicate == "A", method == "Damage_seq", product == "CPD")

# filtering for A.2
pA2_data <- filter(sim_df_org, phase != "async", replicate == "A", 
                   product == "CPD", method == "Damage_seq")

# filtering for B.1
pB1_data <- filter(real_df_org, method != "DNA_seq", phase != "async", 
                   replicate == "A", method == "XR_seq", product == "CPD")

# filtering for B.2
pB2_data <- filter(sim_df_org, phase != "async", method != "DNA_seq", 
                   replicate == "A", product == "CPD", method == "XR_seq")

# for naming of simulated samples
method_labs_sim <- c("Simulated \nDamage-seq", "Simulated \nXR-seq")
names(method_labs_sim) <- c("Damage_seq", "XR_seq")

#### Plot A.1 ####

# create the plot 
p.A.1 <- ggplot(pA1_data, aes(x = windows, y = RPKM)) + 
  geom_vline(xintercept = 0, color = "gray", linetype = "dashed") +
  geom_line(aes(color = sample_strand)) + 
  facet_grid(~product~time_after_exposure~phase~method, 
             labeller = labeller(product = product_labs, 
                                 method = method_labs, 
                                 time_after_exposure = taex_labs,
                                 phase = phase_labs)) + 
  xlab("") + ylab(fr_lab) +
  scale_y_continuous(breaks = c(.1, .2, .3),
                     limits = c(.1, .3)) +
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
  facet_grid(~product~time_after_exposure~phase~method,
             labeller = labeller(product = product_labs,
                                 method = method_labs_sim,
                                 time_after_exposure = taex_labs,
                                 phase = phase_labs)) +
  xlab("") + ylab(fr_lab) +
  scale_y_continuous(breaks = c(.1, .2, .3),
                     limits = c(.1, .3)) +
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
  facet_grid(~product~time_after_exposure~phase~method, 
             labeller = labeller(product = product_labs, 
                                 method = method_labs, 
                                 time_after_exposure = taex_labs,
                                 phase = phase_labs)) + 
  xlab("") + ylab(fr_lab) +
  scale_y_continuous(breaks = c(.1, .2, .3),
                     limits = c(.1, .35)) +
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
  facet_grid(~product~time_after_exposure~phase~method, 
             labeller = labeller(product = product_labs, 
                                 method = method_labs_sim, 
                                 time_after_exposure = taex_labs, 
                                 phase = phase_labs)) + 
  xlab("") + ylab(fr_lab) +
  scale_y_continuous(breaks = c(.1, .2, .3),
                     limits = c(.1, .35)) +
  scale_x_continuous(limits = c(-101, 101), 
                     breaks = c(-101, 0, 101), 
                     labels = c("-10", "0", "+10")) + 
  scale_color_manual(values = strand_colors) + 
  labs(color = "Strands")

# adding and overriding the default plot format
p.B.2 <- p.B.2 + p_format + 
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.border = element_rect(fill = NA))


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


ggsave("~/Desktop/supfig7_new.png", width = 22, height = 18, units = "cm")

