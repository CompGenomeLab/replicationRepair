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

# plus stranded genes csv file for plot A 
gene_plus_csv <- paste("~/Documents/myprojects/replicationRepair/3_output/",
                       "gitignore/1_TextforPlotting/", 
                       "[2020.11.28]plus_genic.csv", 
                       sep = "")

# minus stranded genes csv file for plot A 
gene_minus_csv <- paste("~/Documents/myprojects/replicationRepair/3_output/",
                        "gitignore/1_TextforPlotting/", 
                        "[2020.11.28]minus_genic.csv", 
                        sep = "")

# intergenic regions csv file for plot A 
intergenic_csv <- paste("~/Documents/myprojects/replicationRepair/3_output/",
                        "gitignore/1_TextforPlotting/", 
                        "[2020.11.28]intergenic.csv", 
                        sep = "")

# name of the sample csv file for plot B and C
sample_csv <- paste("~/Documents/myprojects/replicationRepair/3_output/",
                    "gitignore/1_TextforPlotting/", 
                    "[2020.10.13]final_report_inZones_repdomains_intergenic_",
                    "windows_201_100_ready.csv", 
                    sep = "")

# path of the default plot format and functions
sourcePath <- "~/Documents/myprojects/replicationRepair/1_code/r/"


#### Default Plot Format ####

source("4_plot_format.R")


#### Fuctions ####

source("4_functions.R")


#### Main ####

# for plot A
plus_genes_df <- read.csv( gene_plus_csv )
plus_genes_df$region <- "Genic (+)"
minus_genes_df <- read.csv( gene_minus_csv )
minus_genes_df$region <- "Genic (-)"
intergenic_df <- read.csv( intergenic_csv )
intergenic_df$region <- "Intergenic"
pA_data <- rbind(plus_genes_df, minus_genes_df, intergenic_df)

pA_data$region <- factor(
  pA_data$region, levels = c("Genic (+)", "Genic (-)", "Intergenic"))

# for plot B and C
sample_df <- read.csv( sample_csv )
df_rr <- repair_rate( sample_df )
df_rr_org <- window_numbering( df_rr, 4, 101 )
df_rr_org <- domain_name( df_rr_org, 1 )
df_rr_org$dataset <- "Initiation Zones"

df_rr_org$sample_strand <- factor(
  df_rr_org$sample_strand, levels = c("+","-"))

# filtering for B.1
pB1_data <- filter(df_rr_org, phase != "async", replicate == "A", 
                   product == "CPD", time_after_exposure == "12", 
                   phase == "early")

# for plot B.2
pB2_data <- rr_boxplot( pB1_data ) 

# for plot B.3
pB3_data <- rr_boxplot_plus_minus( pB1_data ) 

# filtering for C.1
pC1_data <- filter(df_rr_org, phase != "async", replicate == "A", 
                   product == "CPD", time_after_exposure == "120", 
                   phase == "early")

# for plot C.2
pC2_data <- rr_boxplot( pC1_data ) 

# for plot C.3
pC3_data <- rr_boxplot_plus_minus( pC1_data ) 


#### Plot A ####

# create the plot
p.A <- ggplot(pA_data, aes(x = position, 
                           y = gene_intergenic_percentage)) +
  geom_smooth(aes(color = region), se=FALSE) +
  xlab("Position Relative to Initiation Zones (kb)") + 
  ylab("(%)") +
  scale_x_continuous(limits = c(-10000, 10000), 
                     breaks = c(-10000, 0, 10000), 
                     labels = c("-10", "0", "+10")) +
  scale_color_manual(values = c("#0571b0", "#ca0020", "chartreuse4")) +
  labs(color = "Regions")

# adding and overriding the default plot format
p.A <- p.A + p_format + 
  theme(legend.position = "right")


#### Plot B.1 ####

# create the plot 
p.B.1 <- ggplot(pB1_data, aes(x = windows, y = log2(xr_ds))) +
  geom_vline(xintercept = 0, color = "gray", linetype = "dashed") + 
  geom_line(aes(color = sample_strand)) + 
  facet_grid(~repdomains) +
  xlab(windows_lab) + ylab(fr_xr_ds_lab) +
  scale_x_continuous(limits = c(-101, 101), 
                     breaks = c(-101, 0, 101), 
                     labels = c("-10", "0", "+10")) + 
  scale_color_manual(values = strand_colors) + 
  ylim(-1.5, 1.5) + 
  guides(color = FALSE)

# adding and overriding the default plot format
p.B.1 <- p.B.1 + p_format + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_text(hjust=c(0.3, 0.5, 0.7)))


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
                     limits = c(0, 2.5)) 

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
  xlab("Position Relative to Initiation Zones (kb)\n(Intergenic)") + 
  ylab(fr_xr_ds_lab) +
  scale_x_continuous(limits = c(-101, 101), 
                     breaks = c(-101, 0, 101), 
                     labels = c("-10", "0", "+10")) + 
  scale_color_manual(values = strand_colors) + 
  ylim(-1.5, 1.5) + 
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

p.B.1.2.3 <- p.B.1 + (p.B.2 / p.B.3) 
p.C.1.2.3 <- p.C.1 + (p.C.2 / p.C.3) +
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")
p.A_new <- p.A + grid::textGrob("") +
  plot_layout(guides = "keep")

layout <- "
AAAABBB
CCCCCCD
CCCCCCD
EEEEEEF
EEEEEEF
"


p.A_new + p.B.1.2.3 + 
  grid::textGrob('CPD\n12 min.\nEarly S Phase', 
                                       rot = -90, gp=gpar(fontsize=12), 
                                       y = unit(.55, "npc")) + 
  p.C.1.2.3 + grid::textGrob('CPD\n120 min.\nEarly S Phase', 
                                   rot = -90, gp=gpar(fontsize=12), 
                                   y = unit(.62, "npc")) + 
  plot_layout(design = layout) & 
  theme(plot.tag = element_text(size = 12, face="bold"), 
        plot.title = element_text(hjust = -0.2, vjust = 5, 
                                  size = 12, face="bold"))


ggsave("~/Desktop/supfig11.svg", width = 22, height = 18, units = "cm")

