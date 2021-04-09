#### Packages and Libraries ####

library(stringr)
library(ggpubr)
library(ggplot2)
library(dplyr)
library(reshape2)
library(patchwork)
library(ggthemes)


#### Variables ####

# name of the mutation bed file 
mut_bed <- paste("~/Documents/myprojects/replicationRepair/3_output/",
                 "gitignore/1_TextforPlotting/", 
                 "[2020.10.05.14:33]inzones_repdomains_201_100_",
                 "melanoma_mutations_combined.bed",
                 sep = "")

# name of the TC dipyrimidine content file 
TCcontent_txt <- paste("~/Documents/myprojects/replicationRepair/3_output/",
                       "gitignore/1_TextforPlotting/", 
                       "[2020.10.19]InZones_nuc_content.txt",
                       sep = "")

# path of the default plot format and functions
sourcePath <- "~/Documents/myprojects/replicationRepair/1_code/r/"


#### Default Plot Format ####

source(paste(sourcePath, "4_plot_format.R", sep = ""))


#### Fuctions ####

source(paste(sourcePath, "4_functions.R", sep = ""))


#### Main ####

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

# for plot A.2
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
pA2_data <- rbind(mut_plus_agg, mut_minus_agg)

# for plot A.3
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
pA3_data <- rbind(mut_plus_agg, mut_minus_agg) 


#### Plot A.1 ####

# create the plot
p.A.1 <- ggplot(mut) + 
  geom_line(aes(x = window_number, y = norm, 
                color = strands)) + 
  facet_grid(~repdomains) +
  xlab("Position Relative to Initiation Zones (kb)") + 
  ylab("Normalized Mutation\nCount") +
  geom_vline(xintercept = 0, color = "gray", linetype = "dashed") +
  scale_x_continuous(limits = c(-100, 100), 
                     breaks = c(-100, 0, 100), 
                     labels = c("-10", "0", "+10")) + 
  scale_color_manual(values = strand_colors) + 
  labs(color = "Strands") 

# adding and overriding the default plot format
p.A.1 <- p.A.1 + p_format + 
  theme(legend.position = "bottom",
        panel.border=element_rect(fill = NA),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.text.x = element_text(hjust=c(0.1, 0.5, 0.9))) 


#### Plot A.2 ####

# create the plot
p.A.2 <- ggplot(data = pA2_data, aes(x = Group.2, y = x, 
                                     fill = Group.3)) + 
  geom_bar(stat = "identity", size = 1.5, 
           position=position_dodge2(padding = 0.05)) +
  facet_wrap(~direction) +
  xlab("Replication Domains") + 
  ylab("Mutation Count\n(MC)") +
  scale_fill_manual(values = c("+" = "#0571b0", 
                               "-" = "#ca0020"), 
                    guide = FALSE) +
  scale_y_continuous(breaks = c(.00, .02, .04),
                     limits = c(0, .05)) +
  labs(color = "Strands", fill = "") 

# adding and overriding the default plot format
p.A.2 <- p.A.2 + p_format + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x = element_blank())


#### Plot A.3 ####

# create the plot
p.A.3 <- ggplot(data = pA3_data, aes(x = Group.2, y = x)) + 
  geom_bar(stat = "identity", position=position_dodge()) +
  facet_wrap(~direction) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  xlab("Replication Domains") + ylab(expression(MC[p] - MC[m])) +
  scale_y_continuous(breaks = c(-.003, 0, .003),
                     limits = c(-.005, .005)) 

# adding and overriding the default plot format
p.A.3 <- p.A.3 + p_format + 
  theme(strip.background = element_blank(),
        strip.text.x = element_blank())  


#### Combining Plots with Patchwork ####

layout <- "
AAABB
AAACC
"

(p.A.1 + p.A.2 + p.A.3) +
  plot_layout(design = layout, guides = "collect", tag_level = 'new') & 
  theme(legend.position = 'bottom') 


############### size is problem 
ggsave("~/Desktop/fig6.svg", width = 22, height = 9, units = "cm") 

