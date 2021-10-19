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
sample_csv <- paste("/Users/azgarian/Desktop/final/", 
                    "final_reports_markers_iz_repdomains_uv_m0.5_hela_windows_21_1000_intergenic.txt",
                    sep = "")


#### Default Plot Format ####

source("/Users/azgarian/Documents/myprojects/replicationRepair/workflow/scripts/4plots/4_plot_format.R")


#### Fuctions ####

source("/Users/azgarian/Documents/myprojects/replicationRepair/workflow/scripts/4plots/4_functions.R")


#### Main ####

sample_df <- read.table( sample_csv )

colnames(sample_df) <- c("chromosomes", "start_position", "end_position", "dataset", 
                  "score", "dataset_strand", "counts", "sample_names", 
                  "sample_strand", 
                  "mapped_reads", "RPKM")

df_rr_org <- window_numbering( sample_df, 4, 11 )
df_rr_org <- domain_name( df_rr_org, 1 )
df_rr_org$dataset <- "Initiation Zones"

df_rr_org$sample_strand <- factor(
  df_rr_org$sample_strand, levels = c("+","-"))


# create the plot 
p <- ggplot(df_rr_org, aes(x = windows, y = RPKM)) +
  geom_vline(xintercept = 0, color = "gray", linetype = "dashed") + 
  geom_line() + 
  facet_grid(~sample_names~repdomains) +
  xlab(windows_lab) + ylab("RPKM") +
  scale_x_continuous(limits = c(-11, 11), 
                     breaks = c(-11, 0, 11), 
                     labels = c("-10", "0", "+10")) + 
  scale_color_manual(values = strand_colors) + 
  labs(color = "Strands")
p
# adding and overriding the default plot format
p.B.1 <- p.B.1 + p_format + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_text(hjust=c(0.3, 0.5, 0.7)))


#### Plot B.2 ####

# "Plus Minus" writing in the plot
dat_text <- data.frame(
  label = c("Plus\nMinus", ""),
  direction = c("Left Replicating", "Right Replicating"),
  x = c(2, 2), y = c(0.8, 0.8))

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
  labs(color = "Strands", fill = "") +
  scale_y_continuous(breaks = c(0, 1, 2),
                     limits = c(0, 2)) +
  geom_text( data = dat_text, mapping = aes(x = x, y = y, label = label), 
             angle = 90, colour = "white", size = 2.8 ) 

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
  xlab("") + ylab(expression(RR[p] - RR[m])) +
  scale_y_continuous(breaks = c(-.2, 0, .2),
                     limits = c(-.25, .25)) +
  scale_fill_manual(values = repdomain_colors, guide = FALSE) 

# adding and overriding the default plot format
p.B.3 <- p.B.3 + p_format + 
  theme(strip.background = element_blank(),
        strip.text.x = element_blank())   


#### Plot C.1 ####

# create the plot 
p.C.1 <- ggplot(pC1_data, aes(x = windows, y = log2(xr_ds))) + 
  geom_vline(xintercept = 0, color = "gray", linetype = "dashed") +
  geom_line(aes(color = sample_strand)) + 
  facet_grid(~repdomains) +
  xlab("Position Relative to Initiation Zones (kb)") + 
  ylab(fr_xr_ds_lab) +
  scale_x_continuous(limits = c(-101, 101), 
                     breaks = c(-101, 0, 101), 
                     labels = c("-10", "0", "+10")) + 
  scale_color_manual(values = strand_colors) + 
  labs(color = "Strands")

# adding and overriding the default plot format
p.C.1 <- p.C.1 + p_format + 
  theme(axis.text.x=element_text(hjust=c(0.3, 0.5, 0.7)))


#### Plot C.2 ####

# create the plot 
p.C.2 <- ggplot() + 
  geom_bar(data = pC2_data, aes(x = Group.2, y = x, 
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
p.C.2 <- p.C.2 + p_format + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x = element_blank()) 


#### Plot C.3 ####

# create the plot 
p.C.3 <- ggplot() + 
  geom_bar(data = pC3_data, aes(x = Group.2, y = x), 
           stat = "identity", position=position_dodge()) +
  facet_wrap(~direction) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  xlab("Replication Domains") + 
  ylab(expression(RR[p] - RR[m])) +
  scale_fill_manual(values = repdomain_colors, guide = FALSE) +
  scale_y_continuous(breaks = c(-.2, 0, .2),
                     limits = c(-.25, .25))

# adding and overriding the default plot format
p.C.3 <- p.C.3 + p_format + 
  theme(strip.background = element_blank(),
        strip.text.x = element_blank()) 


#### Combining Plots with Patchwork ####

p.B.2.3 <- (p.B.2 / p.B.3) 
p.C.2.3 <- (p.C.2 / p.C.3)  

layout <- "
AAAAAAA
BBBCCCD
FFFGGGH
"

p.A + p.B.1 + p.B.2.3 + grid::textGrob('CPD\n12 min.\nEarly S Phase', 
                                       rot = -90, gp=gpar(fontsize=12), 
                                       y = unit(.55, "npc")) + 
  p.C.1 + p.C.2.3 + grid::textGrob('CPD\n120 min.\nEarly S Phase', 
                                   rot = -90, gp=gpar(fontsize=12), 
                                   y = unit(.62, "npc")) + 
  plot_layout(design = layout, guides = "collect") & 
  theme(plot.tag = element_text(size = 12, face="bold"),
        legend.position = 'bottom', 
        plot.title = element_text(hjust = -0.2, vjust = 5, 
                                  size = 12, face="bold"))




ggsave("~/Desktop/iz_uv_early_10kb.png", width = 22, height = 18, units = "cm")


