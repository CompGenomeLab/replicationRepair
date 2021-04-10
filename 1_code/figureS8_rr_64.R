#### Packages and Libraries ####

library(stringr)
library(ggpubr)
library(ggplot2)
library(dplyr)
library(reshape2)
library(patchwork)
library(ggthemes)
library(grid)


#### Variables ####

# name of the sample csv file 
sample_csv <- paste("~/Documents/myprojects/replicationRepair/3_output/",
                    "gitignore/1_TextforPlotting/", 
                    "[2020.02.10]final_report_inZones_repdomains_",
                    "windows_201_100_ready.csv", 
                    sep = "")


#### Default Plot Format ####

source("4_plot_format.R")


#### Fuctions ####

source("4_functions.R")


#### Main ####

sample_df <- read.csv( sample_csv )
df_rr <- repair_rate( sample_df )
df_rr_org <- window_numbering( df_rr, 4, 101 )
df_rr_org <- domain_name( df_rr_org, 1 )
df_rr_org$dataset <- "Initiation Zones"

df_rr_org$sample_strand <- factor(
  df_rr_org$sample_strand, levels = c("+","-"))

# filtering for A.1
pA1_data <- filter(df_rr_org, phase != "async", replicate == "A", 
                   product == "64_PP", time_after_exposure == "12", 
                   phase == "early")

# for plot A.2
pA2_data <- rr_boxplot( pA1_data ) 

# for plot A.3
pA3_data <- rr_boxplot_plus_minus( pA1_data ) 

# filtering for B.1
pB1_data <- filter(df_rr_org, phase != "async", replicate == "A", 
                   product == "64_PP", time_after_exposure == "12", 
                   phase == "late")

# for plot B.2
pB2_data <- rr_boxplot( pB1_data ) 

# for plot B.3
pB3_data <- rr_boxplot_plus_minus( pB1_data ) 


#### Plot A.1 ####

# create the plot 
p.A.1 <- ggplot(pA1_data, aes(x = windows, y = log2(xr_ds))) +
  geom_vline(xintercept = 0, color = "gray", linetype = "dashed") + 
  geom_line(aes(color = sample_strand)) + 
  facet_grid(~repdomains) +
  xlab(windows_lab) + ylab(fr_xr_ds_lab) +
  scale_x_continuous(limits = c(-101, 101), 
                     breaks = c(-101, 0, 101), 
                     labels = c("-10", "0", "+10")) + 
  scale_color_manual(values = strand_colors) + 
  ylim(-1.5, 1.5) + 
  labs(color = "Strands")

# adding and overriding the default plot format
p.A.1 <- p.A.1 + p_format + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_text(hjust=c(0.3, 0.5, 0.7)))


#### Plot A.2 ####

# create the plot 
p.A.2 <- ggplot() + 
  geom_bar(data = pA2_data, aes(x = Group.2, y = x, 
                                fill = Group.3), 
           stat = "identity", size = 1.5, 
           position=position_dodge2(padding = 0.05)) +
  facet_wrap(~direction) +
  xlab("Replication Domains") + 
  ylab("Repair\nRate (RR)") +
  scale_y_continuous(breaks = c(0, 1, 2),
                     limits = c(0, 2)) +
  scale_fill_manual(values = c("+" = "#0571b0", 
                               "-" = "#ca0020"), 
                    guide = FALSE) +
  labs(color = "Strands", fill = "") 

# adding and overriding the default plot format
p.A.2 <- p.A.2 + p_format + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x = element_blank()) 


#### Plot A.3 ####

# create the plot 
p.A.3 <- ggplot() + 
  geom_bar(data = pA3_data, aes(x = Group.2, y = x), 
           stat = "identity", position=position_dodge()) +
  facet_wrap(~direction) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  xlab("") + ylab(expression(RR[p] - RR[m])) +
  scale_y_continuous(breaks = c(-.2, 0, .2),
                     limits = c(-.25, .25)) +
  scale_fill_manual(values = repdomain_colors, guide = FALSE) 

# adding and overriding the default plot format
p.A.3 <- p.A.3 + p_format + 
  theme(strip.background = element_blank(),
        strip.text.x = element_blank())   


#### Plot B.1 ####

# create the plot 
p.B.1 <- ggplot(pB1_data, aes(x = windows, y = log2(xr_ds))) + 
  geom_vline(xintercept = 0, color = "gray", linetype = "dashed") +
  geom_line(aes(color = sample_strand)) + 
  facet_grid(~repdomains) +
  xlab("Position Relative to Initiation Zones (kb)") + 
  ylab(fr_xr_ds_lab) +
  scale_x_continuous(limits = c(-101, 101), 
                     breaks = c(-101, 0, 101), 
                     labels = c("-10", "0", "+10")) + 
  scale_color_manual(values = strand_colors) + 
  ylim(-1.5, 1.5) + 
  labs(color = "Strands")

# adding and overriding the default plot format
p.B.1 <- p.B.1 + p_format + 
  theme(axis.text.x=element_text(hjust=c(0.3, 0.5, 0.7)))


#### Plot B.2 ####

# create the plot 
p.B.2 <- ggplot() + 
  geom_bar(data = pB2_data, aes(x = Group.2, y = x, 
                                fill = Group.3), 
           stat = "identity", size = 1.5, 
           position=position_dodge2(padding = 0.05)) +
  facet_wrap(~direction) +
  xlab("Replication Domains") + 
  ylab("Repair\nRate (RR)") +
  scale_fill_manual(values = c("+" = "#0571b0", 
                               "-" = "#ca0020"), 
                    guide = FALSE) +
  scale_y_continuous(breaks = c(0, 1, 2),
                     limits = c(0, 2)) +
  labs(color = "Strands", fill = "") 

# adding and overriding the default plot format
p.B.2 <- p.B.2 + p_format + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x = element_blank()) 


#### Plot B.3 ####

# create the plot 
p.B.3 <- ggplot() + 
  geom_bar(data = pB3_data, aes(x = Group.2, y = x), 
           stat = "identity", position=position_dodge()) +
  facet_wrap(~direction) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  xlab("Replication Domains") + 
  ylab(expression(RR[p] - RR[m])) +
  scale_y_continuous(breaks = c(-.2, 0, .2),
                     limits = c(-.25, .25)) +
  scale_fill_manual(values = repdomain_colors, guide = FALSE) 

# adding and overriding the default plot format
p.B.3 <- p.B.3 + p_format + 
  theme(strip.background = element_blank(),
        strip.text.x = element_blank()) 


#### Combining Plots with Patchwork ####

p.A.2.3 <- (p.A.2 / p.A.3)  
p.B.2.3 <- (p.B.2 / p.B.3) 


layout <- "
BBBCCCD
FFFGGGH
"

p.A.1 + p.A.2.3 + grid::textGrob(('(6-4)PP\n12 min.\nEarly S Phase'), 
                                 rot = -90, gp=gpar(fontsize=12), 
                                 y = unit(.55, "npc")) + 
  p.B.1 + p.B.2.3 + grid::textGrob(('(6-4)PP\n12 min.\nLate S Phase'), 
                                   rot = -90, gp=gpar(fontsize=12), 
                                   y = unit(.62, "npc")) + 
  plot_layout(design = layout, guides = "collect") & 
  theme(plot.tag = element_text(size = 12, face="bold"),
        legend.position = 'bottom', 
        plot.title = element_text(hjust = -0.2, vjust = 5, 
                                  size = 12, face="bold"))




ggsave("~/Desktop/supfig8.png", width = 22, height = 18, units = "cm")
