#### Packages and Libraries ####

library(stringr)
library(ggpubr)
library(ggplot2)
library(dplyr)
library(reshape2)
library(patchwork)
library(ggthemes)


#### Variables ####

# name of the noDrivers mutation bed file 
mut_noDriver_bed <- paste("~/Documents/myprojects/replicationRepair/3_output/",
                     "gitignore/1_TextforPlotting/", 
                     "[2020.10.05.14:34]inzones_repdomains_201_100_melanoma_",
                     "mutations_combined_noDrivers.bed",
                     sep = "")

# name of the noDrivers TC dipyrimidine content file 
TCcontent_noDriver_txt <- paste("~/Documents/myprojects/replicationRepair/",
                                "3_output/", "gitignore/1_TextforPlotting/", 
                                "[2020.10.19]InZones_noDrivers_nuc_content.txt",
                                sep = "")

# name of the intergenic mutation bed file 
mut_intergenic_bed <- paste("~/Documents/myprojects/replicationRepair/",
                            "3_output/gitignore/1_TextforPlotting/", 
                            "[2020.10.13.19:18]inzones_repdomains_201_100_",
                            "melanoma_mutations_combined_intergenic.bed",
                            sep = "")

# name of the intergenic TC dipyrimidine content file 
TCcontent_intergenic_txt <- paste("~/Documents/myprojects/replicationRepair/",
                                  "3_output/gitignore/1_TextforPlotting/", 
                                  "[2020.10.19]InZones_intergenic_nuc_",
                                  "content.txt", sep = "")

p1_mut_name <- "Normalized Mutation\nCount"
p2_mut_name <- "Mutation\nCount (MC)"  


# path of the default plot format and functions
sourcePath <- "~/Documents/myprojects/replicationRepair/1_code/r/"


#### Default Plot Format ####

source(paste(sourcePath, "4_plot_format.R", sep = ""))


#### Fuctions ####

source(paste(sourcePath, "4_functions.R", sep = ""))

TC_content <- function (TCcontent_txt) {
  # TC content df arrangement
  TCcontent <- read.table( TCcontent_txt, header = TRUE)
  TCcontent <- TCcontent[,c(1,6,7)]
  names(TCcontent) <- c("name","+", "-")
  TCcontent_plus <- TCcontent[,c(1,2)]
  TCcontent_minus <- TCcontent[,c(1,3)]
  TCcontent_plus <- melt(TCcontent_plus, measure.vars = "+")
  TCcontent_minus <- melt(TCcontent_minus, measure.vars = "-")
  TCcontent <- rbind(TCcontent_plus, TCcontent_minus)
  names(TCcontent) <- c("name","strands", "TC_content")
  return(TCcontent)
}

mut_arrange <- function (mut_bed, TCcontent) {
  # mutations
  mut_df <- read.table( mut_bed )
  names(mut_df) <- c("chr", "start_pos", "end_pos", "name", "score", "strands", 
                     "count")
  mut <- merge(mut_df, TCcontent, by.x=c("name", "strands"), 
               by.y=c("name", "strands"))
  Windows = data.frame(str_split_fixed(mut$name, "_", -1))
  mut$window_number = as.numeric(levels(Windows[ , ncol(Windows)]
  ))[Windows[ , ncol(Windows)]] -101
  mut$repdomains <- Windows[,3]
  mut$name <- "Initiation_Zones"
  mut$norm <- mut$count / mut$TC_content
  
  mut$repdomains <- factor(mut$repdomains, levels = c("UTZ", "ERD", "DTZ", "LRD"))
  mut$strands <- factor(
    mut$strands, levels = c("+","-"))
  return(mut)
}

plot2_arrange <- function (mut) {
  mut_plus <- filter(mut, window_number > 0)
  mut_plus_agg <- aggregate(x = mut_plus$norm, by = 
                              list(mut_plus$name, mut_plus$repdomains, 
                                   mut_plus$strands), FUN = "mean")
  mut_plus_agg$direction <- "Right Replicating"
  mut_minus <- filter(mut, window_number < 0)
  mut_minus_agg <- aggregate(x = mut_minus$norm, by = 
                               list(mut_minus$name, mut_minus$repdomains, 
                                    mut_minus$strands), FUN = "mean")
  mut_minus_agg$direction <- "Left Replicating"
  p2_data <- rbind(mut_plus_agg, mut_minus_agg)
  return(p2_data)
}

plot3_arrange <- function (mut) {
  mut_casted <- dcast(mut, name + repdomains + window_number ~ strands, 
                      value.var = "norm")
  mut_casted$plus_min <- mut_casted$"+" - mut_casted$"-" 
  mut_plus <- filter(mut_casted, window_number > 0)
  mut_plus_agg <- aggregate(x = mut_plus$plus_min, by = 
                              list(mut_plus$name, mut_plus$repdomains), 
                            FUN = "mean")
  mut_plus_agg$direction <- "Right Replicating"
  mut_minus <- filter(mut_casted, window_number < 0)
  mut_minus_agg <- aggregate(x = mut_minus$plus_min, by = 
                               list(mut_minus$name, mut_minus$repdomains), 
                             FUN = "mean")
  mut_minus_agg$direction <- "Left Replicating"
  p3_data <- rbind(mut_plus_agg, mut_minus_agg) 
  return(p3_data)
}


#### Main ####

# for plot B
pB_tc <- TC_content(TCcontent_noDriver_txt)
pB1_data <- mut_arrange(mut_noDriver_bed, pB_tc)
pB2_data <- plot2_arrange(pB1_data)
pB3_data <- plot3_arrange(pB1_data)

# for plot C
pC_tc <- TC_content(TCcontent_intergenic_txt)
pC1_data <- mut_arrange(mut_intergenic_bed, pC_tc)
pC2_data <- plot2_arrange(pC1_data)
pC3_data <- plot3_arrange(pC1_data)


#### Plot B.1 ####

# create the plot
p.B.1 <- ggplot(pB1_data) + 
  geom_line(aes(x = window_number, y = norm, 
                color = strands)) + 
  facet_grid(~repdomains) +
  xlab("Position Relative to Initiation Zones (kb)\n(without Drivers)") + 
  ylab(p1_mut_name) +
  geom_vline(xintercept = 0, color = "gray", linetype = "dashed") +
  scale_y_continuous(breaks = c(0, .02, .04),
                     limits = c(0, .055)) +
  scale_x_continuous(limits = c(-100, 100), 
                     breaks = c(-100, 0, 100), 
                     labels = c("-10", "0", "+10")) + 
  scale_color_manual(values = strand_colors) + 
  labs(color = "Strands") 

# adding and overriding the default plot format
p.B.1 <- p.B.1 + p_format + 
  theme(legend.position = "bottom",
        panel.border=element_rect(fill = NA),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.text.x = element_text(hjust=c(0.1, 0.5, 0.9))) 


#### Plot B.2 ####

# create the plot
p.B.2 <- ggplot(data = pB2_data, aes(x = Group.2, y = x, 
                                     fill = Group.3)) + 
  geom_bar(stat = "identity", size = 1.5, 
           position=position_dodge2(padding = 0.05)) +
  facet_wrap(~direction) +
  xlab("") + 
  ylab(p2_mut_name) +
  scale_y_continuous(breaks = c(0, .02, .04),
                     limits = c(0, .05)) +
  scale_fill_manual(values = c("+" = "#0571b0", 
                               "-" = "#ca0020"), 
                    guide = FALSE) +
  labs(color = "Strands", fill = "") 

# adding and overriding the default plot format
p.B.2 <- p.B.2 + p_format + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x = element_blank()) 


#### Plot B.3 ####

# create the plot
p.B.3 <- ggplot(data = pB3_data, aes(x = Group.2, y = x)) + 
  geom_bar(stat = "identity", position=position_dodge()) +
  facet_wrap(~direction) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  xlab("") + ylab(expression(MC[p] - MC[m])) +
  scale_y_continuous(breaks = c(-.004, 0, .004),
                     limits = c(-.005, .005))

# adding and overriding the default plot format
p.B.3 <- p.B.3 + p_format + 
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.title.x = element_blank()) 

#### Plot C.1 ####

# create the plot
p.C.1 <- ggplot(pC1_data) + 
  geom_line(aes(x = window_number, y = norm, 
                color = strands)) + 
  facet_grid(~repdomains) +
  xlab("Position Relative to Initiation Zones (kb)\n(intergenic)") + 
  ylab(p1_mut_name) +
  geom_vline(xintercept = 0, color = "gray", linetype = "dashed") +
  scale_y_continuous(breaks = c(0, .02, .04),
                     limits = c(0, .055)) +
  scale_x_continuous(limits = c(-100, 100), 
                     breaks = c(-100, 0, 100), 
                     labels = c("-10", "0", "+10")) + 
  scale_color_manual(values = strand_colors) + 
  labs(color = "Strands") 

# adding and overriding the default plot format
p.C.1 <- p.C.1 + p_format + 
  theme(legend.position = "bottom",
        panel.border=element_rect(fill = NA),
        axis.title.x = element_text(vjust = -1.5),
        axis.text.x = element_text(hjust=c(0.1, 0.5, 0.9)))


#### Plot C.2 ####

# create the plot
p.C.2 <- ggplot(data = pC2_data, aes(x = Group.2, y = x, 
                                     fill = Group.3)) + 
  geom_bar(stat = "identity", size = 1.5, 
           position=position_dodge2(padding = 0.05)) +
  facet_wrap(~direction) +
  xlab("Replication Domains") + 
  ylab(p2_mut_name) +
  scale_y_continuous(breaks = c(0, .02, .04),
                     limits = c(0, .05)) +
  scale_fill_manual(values = c("+" = "#0571b0", 
                               "-" = "#ca0020"), 
                    guide = FALSE) +
  labs(color = "Strands", fill = "") 

# adding and overriding the default plot format
p.C.2 <- p.C.2 + p_format + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x = element_blank())


#### Plot C.3 ####

# create the plot
p.C.3 <- ggplot(data = pC3_data, aes(x = Group.2, y = x)) + 
  geom_bar(stat = "identity", position=position_dodge()) +
  facet_wrap(~direction) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  xlab("Replication Domains") + ylab(expression(MC[p] - MC[m])) +
  scale_y_continuous(breaks = c(-.004, 0, .004),
                     limits = c(-.005, .005))

# adding and overriding the default plot format
p.C.3 <- p.C.3 + p_format + 
  theme(strip.background = element_blank(),
        axis.title.x = element_text(vjust = -1.5),
        strip.text.x = element_blank()) 


#### Combining Plots with Patchwork ####

layout <- "
DDDEE
DDDFF
GGGHH
GGGII
"
options(scipen=999)

(p.B.1 + p.B.2 + p.B.3 + p.C.1 + p.C.2 + p.C.3) +
  plot_layout(design = layout, guides = "collect") & 
  plot_annotation(tag_levels = 'supp_muts') &
  theme(plot.tag = element_text(size = 12, face="bold"),
        legend.position = 'bottom', 
        plot.title = element_text(hjust = -0.2, vjust = 5, 
                                  size = 12, face="bold"))

ggsave("~/Desktop/supfig12.svg", width = 22, height = 18, units = "cm") 

