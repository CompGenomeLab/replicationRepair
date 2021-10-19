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
sample_sim_csv <- paste("~/Documents/myprojects/replicationRepair/results/final/", 
                    "final_reports_sim_hg19_iz_repdomains_m0.5_hela",
                    "_windows_201_100_intergenic_ready.csv", 
                    sep = "")

sample_csv <- paste("~/Documents/myprojects/replicationRepair/results/final/", 
                        "final_reports_hg19_iz_repdomains_m0.5_hela",
                        "_windows_201_100_intergenic_ready.csv", 
                        sep = "")

#### Default Plot Format ####

source("~/Documents/myprojects/replicationRepair/workflow/scripts/4plots/4_plot_format.R")


#### Fuctions ####

source("~/Documents/myprojects/replicationRepair/workflow/scripts/4plots/4_functions.R")

real_sim <- function ( df, sim_rep = "A" ){
  sim <- filter(df, type == "sim" & replicate == sim_rep)
  real <- filter(df, type == "real")
  real_list <- split(real, real$replicate)
  df_real_sim <- real_list[[1]][0, ]
  
  for ( real_rep in 1:length(real_list) ){
    temp <- real_list[[real_rep]]
    temp <- rbind(temp, sim) 
    temp <- dcast(temp, chromosomes + start_position + end_position + 
                    dataset + score + dataset_strand + product + phase + windows +
                    time_after_exposure + repdomains + sample_strand ~ type, 
                  value.var = "xr_ds")
    temp <- cbind(temp, real_list[[real_rep]]["replicate"])
    df_real_sim <- rbind(df_real_sim, temp)
  }
  
  df_real_sim$real_sim <- df_real_sim$real / df_real_sim$sim
  df_real_sim <- select(df_real_sim, -c("real", "sim"))
  
  return(df_real_sim)
}

rr_boxplot <- function (df){
  mut_plus <- filter(df, windows > 0)
  mut_plus_agg <- aggregate(x = mut_plus$real_sim, by = 
                              list(mut_plus$dataset, mut_plus$repdomains, 
                                   mut_plus$sample_strand), FUN = "mean")
  mut_plus_agg$direction <- "Right Replicating"
  mut_minus <- filter(df, windows < 0)
  mut_minus_agg <- aggregate(x = mut_minus$real_sim, by = 
                               list(mut_minus$dataset, mut_minus$repdomains, 
                                    mut_minus$sample_strand), FUN = "mean")
  mut_minus_agg$direction <- "Left Replicating"
  mut_agg <- rbind(mut_plus_agg, mut_minus_agg)
  
  return(mut_agg)
}

rr_boxplot_plus_minus <- function (df){
  mut_casted <- dcast(df, dataset + repdomains + windows ~ sample_strand, 
                      value.var = "real_sim")
  mut_casted$plus_min <- mut_casted$"+" - mut_casted$"-" 
  mut_plus <- filter(mut_casted, windows > 0)
  mut_plus_agg <- aggregate(x = mut_plus$plus_min, by = 
                              list(mut_plus$dataset, mut_plus$repdomains), 
                            FUN = "mean")
  mut_plus_agg$direction <- "Right Replicating"
  mut_minus <- filter(mut_casted, windows < 0)
  mut_minus_agg <- aggregate(x = mut_minus$plus_min, by = 
                               list(mut_minus$dataset, mut_minus$repdomains), 
                             FUN = "mean")
  mut_minus_agg$direction <- "Left Replicating"
  mut_agg <- rbind(mut_plus_agg, mut_minus_agg)
  
  return(mut_agg)
}


#### Main ####

sample_df <- read.csv( sample_csv )
df_rr <- repair_rate( sample_df )
df_rr_org <- window_numbering( df_rr, 4, 101 )
df_rr_org <- domain_name( df_rr_org, 1 )
df_rr_org$dataset <- "Initiation Zones"

df_rr_org$sample_strand <- factor(
  df_rr_org$sample_strand, levels = c("+","-"))

sample_df_sim <- read.csv( sample_sim_csv )
df_rr_sim <- repair_rate( sample_df_sim )
df_rr_org_sim <- window_numbering( df_rr_sim, 4, 101 )
df_rr_org_sim <- domain_name( df_rr_org_sim, 1 )
df_rr_org_sim$dataset <- "Initiation Zones"

df_rr_org_sim$sample_strand <- factor(
  df_rr_org_sim$sample_strand, levels = c("+","-"))

df_rr_org$type <- "real"
df_rr_org_sim$type <- "sim"

df_rr_comb <- rbind(df_rr_org, df_rr_org_sim)

df_rr_rs <- real_sim(df_rr_comb)


# filtering for B.1
pB1_data <- filter(df_rr_rs, phase != "async", replicate == "A", 
                   product == "CPD", time_after_exposure == "12", 
                   phase == "early")

# for plot B.2
pB2_data <- rr_boxplot( pB1_data ) 

# for plot B.3
pB3_data <- rr_boxplot_plus_minus( pB1_data ) 

# filtering for C.1
pC1_data <- filter(df_rr_rs, phase != "async", replicate == "A", 
                   product == "CPD", time_after_exposure == "120", 
                   phase == "early")

# for plot C.2
pC2_data <- rr_boxplot( pC1_data ) 

# for plot C.3
pC3_data <- rr_boxplot_plus_minus( pC1_data ) 


#### Plot A ####

# plot A will be a drawing
# we are creating an empty text for that part
p.A <- wrap_elements(grid::textGrob(''))


#### Plot B.1 ####

# create the plot 
p.B.1 <- ggplot(pB1_data, aes(x = windows, y = log2(real_sim))) +
  geom_vline(xintercept = 0, color = "gray", linetype = "dashed") + 
  geom_line(aes(color = sample_strand)) + 
  facet_grid(~repdomains) +
  xlab(windows_lab) + ylab("normalized Repair Rates (log2)") +
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
  ylab("n. Repair\nRate (nRR)") +
  scale_fill_manual(values = c("+" = "#0571b0", 
                               "-" = "#ca0020"), 
                    guide = "none") +
  labs(color = "Strands", fill = "") +
  #scale_y_continuous(breaks = c(0, 1, 2),
  #                   limits = c(0, 2)) +
  geom_text( data = dat_text, mapping = aes(x = x, y = y, label = label), 
             angle = 90, colour = "white", size = 2.8 ) 

# adding and overriding the default plot format
p.B.2 <- p.B.2 + p_format + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x = element_blank()) 

pB1_data$direction[pB1_data$windows<0] <- "Left Replicating" 
pB1_data$direction[pB1_data$windows>0] <- "Right Replicating" 

pB1_data <- filter(pB1_data, direction != "NA")
pB1_data$log2val <- log2(pB1_data$real_sim)

pB1_boxplot <- pB1_data[,c("repdomains", "log2val", "sample_strand", "direction", "windows")]

p.B.2.v2 <- ggplot(data=pB1_boxplot, aes(x=sample_strand, y=log2val, fill=sample_strand)) +
  geom_boxplot() +
  geom_line(aes(group=windows), colour="grey", size=0.4) +
  facet_wrap(~direction~repdomains) +
  xlab("Replication Domains") + 
  ylab("n. Repair\nRate (nRR)") +
  scale_fill_manual(values = c("+" = "#0571b0", 
                               "-" = "#ca0020"), 
                    guide = "none") +
  #scale_y_continuous(breaks = c(0, 1, 2),
  #                   limits = c(0, 2)) +
  labs(color = "Strands", fill = "") 

# adding and overriding the default plot format
p.B.2.v2 <- p.B.2.v2 + p_format + 
  stat_compare_means(label = "p.format", paired = TRUE) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x = element_blank()) 

#p.B.2.v2

p.B.2.v3 <- ggplot(data=pB1_boxplot, aes(x=repdomains, y=log2val, fill=sample_strand)) +
  geom_boxplot() +
  facet_wrap(~direction) +
  xlab("Replication Domains") + 
  ylab("n. Repair\nRate (nRR)") +
  scale_fill_manual(values = c("+" = "#0571b0", 
                               "-" = "#ca0020"), 
                    guide = "none") +
  #scale_y_continuous(breaks = c(0, 1, 2),
  #                   limits = c(0, 2)) +
  labs(color = "Strands", fill = "") 

# adding and overriding the default plot format
p.B.2.v3 <- p.B.2.v3 + p_format + 
  stat_compare_means(label = "p.signif", label.y= 0,  paired = TRUE) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x = element_blank()) 

#p.B.2.v3

#### Plot B.3 ####

# create the plot 
p.B.3 <- ggplot() + 
  geom_bar(data = pB3_data, aes(x = Group.2, y = x), 
           stat = "identity", position=position_dodge()) +
  facet_wrap(~direction) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  xlab("") + ylab(expression(nRR[p] - nRR[m])) +
  scale_y_continuous(breaks = c(-.2, 0, .2),
                     limits = c(-.25, .25)) +
  scale_fill_manual(values = repdomain_colors, guide = "none") 

# adding and overriding the default plot format
p.B.3 <- p.B.3 + p_format + 
  theme(strip.background = element_blank(),
        strip.text.x = element_blank())   


#### Plot C.1 ####

# create the plot 
p.C.1 <- ggplot(pC1_data, aes(x = windows, y = log2(real_sim))) + 
  geom_vline(xintercept = 0, color = "gray", linetype = "dashed") +
  geom_line(aes(color = sample_strand)) + 
  facet_grid(~repdomains) +
  xlab("Position Relative to Initiation Zones (kb)") + 
  ylab("normalized Repair Rates (log2)") +
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
  ylab("n. Repair\nRate (nRR)") +
  scale_fill_manual(values = c("+" = "#0571b0", 
                               "-" = "#ca0020"), 
                    guide = "none") +
  #scale_y_continuous(breaks = c(0, 1, 2),
  #                   limits = c(0, 2)) +
  labs(color = "Strands", fill = "") 

# adding and overriding the default plot format
p.C.2 <- p.C.2 + p_format + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x = element_blank()) 

pC1_data$direction[pC1_data$windows<0] <- "Left Replicating" 
pC1_data$direction[pC1_data$windows>0] <- "Right Replicating" 

pC1_data <- filter(pC1_data, direction != "NA")
pC1_data$log2val <- log2(pC1_data$real_sim)

pC1_boxplot <- pC1_data[,c("repdomains", "log2val", "sample_strand", "direction", "windows")]

p.C.2.v2 <- ggplot(data=pC1_boxplot, aes(x=sample_strand, y=log2val, fill=sample_strand)) +
  geom_boxplot() +
  geom_line(aes(group=windows), colour="grey", size=0.4) +
  facet_wrap(~direction~repdomains) +
  xlab("Replication Domains") + 
  ylab("n. Repair\nRate (nRR)") +
  scale_fill_manual(values = c("+" = "#0571b0", 
                               "-" = "#ca0020"), 
                    guide = "none") +
  #scale_y_continuous(breaks = c(0, 1, 2),
  #                   limits = c(0, 2)) +
  labs(color = "Strands", fill = "") 

# adding and overriding the default plot format
p.C.2.v2 <- p.C.2.v2 + p_format + 
  stat_compare_means(label = "p.format", paired = TRUE) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x = element_blank()) 

#p.C.2.v2

p.C.2.v3 <- ggplot(data=pC1_boxplot, aes(x=repdomains, y=log2val, fill=sample_strand)) +
  geom_boxplot() +
  facet_wrap(~direction) +
  xlab("Replication Domains") + 
  ylab("n. Repair\nRate (nRR)") +
  scale_fill_manual(values = c("+" = "#0571b0", 
                               "-" = "#ca0020"), 
                    guide = "none") +
  #scale_y_continuous(breaks = c(0, 1, 2),
  #                   limits = c(0, 2)) +
  labs(color = "Strands", fill = "") 

# adding and overriding the default plot format
p.C.2.v3 <- p.C.2.v3 + p_format + 
  stat_compare_means(label = "p.signif", label.y= 0,  paired = TRUE) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x = element_blank()) 

#p.C.2.v3

#### Plot C.3 ####

# create the plot 
p.C.3 <- ggplot() + 
  geom_bar(data = pC3_data, aes(x = Group.2, y = x), 
           stat = "identity", position=position_dodge()) +
  facet_wrap(~direction) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  xlab("Replication Domains") + 
  ylab(expression(nRR[p] - nRR[m])) +
  scale_fill_manual(values = repdomain_colors, guide = "none") +
  scale_y_continuous(breaks = c(-.2, 0, .2),
                     limits = c(-.25, .25))

# adding and overriding the default plot format
p.C.3 <- p.C.3 + p_format + 
  theme(strip.background = element_blank(),
        strip.text.x = element_blank()) 


#### Combining Plots with Patchwork ####

layout <- "
A
A
B
"

p.B.2.3 <- (p.B.2.v3 / p.B.3) + plot_layout(design = layout) 
p.C.2.3 <- (p.C.2.v3 / p.C.3) + plot_layout(design = layout) 

layout2 <- "
AAAAAAA
BBBCCCD
FFFGGGH
"
layout3 <- "
BBBCCCD
FFFGGGH
"
p.B.1 + p.B.2.3 + grid::textGrob('CPD\n12 min.\nEarly S Phase', 
                                 rot = -90, gp=gpar(fontsize=12), 
                                 y = unit(.55, "npc")) + 
  p.C.1 + p.C.2.3 + grid::textGrob('CPD\n120 min.\nEarly S Phase', 
                                   rot = -90, gp=gpar(fontsize=12), 
                                   y = unit(.62, "npc")) + 
  plot_layout(design = layout3, guides = "collect") & 
  theme(plot.tag = element_text(size = 12, face="bold"),
        legend.position = 'bottom', 
        plot.title = element_text(hjust = -0.2, vjust = 5, 
                                  size = 12, face="bold"))




ggsave("~/Desktop/iz_real_sim_early_10kb_deneme.png", width = 22, height = 18, units = "cm")
