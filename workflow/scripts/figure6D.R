#### Packages and Libraries ####

library(argparser)
library(stringr)
library(ggpubr)
library(ggplot2)
library(dplyr)
library(reshape2)
library(patchwork)

######## Arguments ##########
p <- arg_parser("producing the figure 6D")
p <- add_argument(p, "--noverlap", help="mutation counts on non-overlapping HeLa initiation zones")
p <- add_argument(p, "--overlap", help="mutation counts on overlapping HeLa initiation zones")
p <- add_argument(p, "--noverlap_TC", help="TC counts on non-overlapping HeLa initiation zones")
p <- add_argument(p, "--overlap_TC", help="TC counts on overlapping HeLa initiation zones")
p <- add_argument(p, "--fig6D", help="figure output")

# Parse the command line arguments
argv <- parse_args(p)

hela2gmimr_bed <- argv$overlap

helanover_bed <- argv$noverlap
  
hela2gmimr_nuc <- argv$overlap_TC
  
helanover_nuc <- argv$noverlap_TC

#hela2gmimr_bed <- "/home/azgarian/Documents/myprojects/replicationRepair/3_output/gitignore/1_TextforPlotting/[2021.03.30.12/44]iz_hela_to_gm_imr_repdomains_201_100_melanoma_mutations_combined_intergenic.bed"
#helanover_bed <- "~/Documents/myprojects/replicationRepair/results/final/melanoma_target_mut_comb_hela_no_overlap_repdomains_windows_201_100_intergenic_combined.txt"
#hela2gmimr_nuc <- "/home/azgarian/Desktop/repairRep_revision/inZones/IZ.hela_to_gm_imr_repdomains_windows_201_100_counts.txt"
#helanover_nuc <- "~/Documents/myprojects/replicationRepair/results/final/hela_no_overlap_repdomains_windows_201_100_counts.txt"

#### Default Plot Format ####

source("workflow/scripts/plot_format.R")


#### Functions ####

source("workflow/scripts/functions.R")


#### Main ####

# TC content df arrangement
TCcontent <- read.table( hela2gmimr_nuc, header = TRUE)

TCcontent$`+` <- TCcontent$tc + TCcontent$cc
TCcontent$`-` <- TCcontent$ga + TCcontent$gg

TCcontent <- TCcontent[,c(1,6,7)]
names(TCcontent) <- c("name","+", "-")
TCcontent_plus <- TCcontent[,c(1,2)]
TCcontent_minus <- TCcontent[,c(1,3)]
TCcontent_plus <- melt(TCcontent_plus, measure.vars = "+")
TCcontent_minus <- melt(TCcontent_minus, measure.vars = "-")
TCcontent <- rbind(TCcontent_plus, TCcontent_minus)
names(TCcontent) <- c("name","strands", "TC_content")

# mutations
mut_df <- read.table( hela2gmimr_bed )

mut_df <- mut_df[,c(1,2,3,4,5,9,11)]

names(mut_df) <- c("chr", "start_pos", "end_pos", "name", "score", "strands", 
                   "count")
mut <- merge(mut_df, TCcontent, by.x=c("name", "strands"), 
             by.y=c("name", "strands"))
window_number = data.frame(str_split_fixed(mut$name, "_", -1))
mut$window_number <- as.numeric(
  as.character(str_split_fixed(mut$name, "_",3)[,3])) - 101
mut$repdomains <- window_number[,2]
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
  xlab("Position Relative to Initiation Zones (kb)\nIntergenic") + 
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

mut$direction[mut$window_number<0] <- "Left Replicating" 
mut$direction[mut$window_number>0] <- "Right Replicating" 

mut <- filter(mut, direction != "NA")
mut$log2val <- log2(mut$norm)

pA1_boxplot <- mut[,c("repdomains", "log2val", "strands", "direction", "window_number")]

p.A.2 <- ggplot(data=pA1_boxplot, aes(x=repdomains, y=log2val, fill=strands)) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~direction) +
  xlab("Replication Domains") + 
  ylab("n. Repair\nRate (RR)") +
  scale_fill_manual(values = c("+" = "#0571b0", 
                               "-" = "#ca0020"), 
                    guide = "none") +
  #scale_y_continuous(breaks = c(0, 1, 2),
  #                   limits = c(0, 2)) +
  labs(color = "Strands", fill = "") 

# adding and overriding the default plot format
p.A.2 <- p.A.2 + p_format + 
  stat_compare_means(label = "p.signif",  paired = TRUE) + 
  theme(strip.background = element_blank(),
        strip.text.x = element_blank()) 

#### Plot A.3 ####

# create the plot
p.A.3 <- ggplot(data = pA3_data, aes(x = Group.2, y = x)) + 
  geom_bar(stat = "identity", position=position_dodge()) +
  facet_wrap(~direction) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  xlab("Replication Domains") + ylab(expression(MC[p] - MC[m])) +
  scale_y_continuous(breaks = c(-.003, 0, .003),
                     limits = c(-.004, .004)) 

# adding and overriding the default plot format
p.A.3 <- p.A.3 + p_format + 
  theme(strip.background = element_blank(),
        strip.text.x = element_blank())  


#### Combining Plots with Patchwork ####

layout <- "
AABB
"

(p.A.1 + p.A.2) +
  plot_layout(design = layout, guides = "collect", tag_level = 'new') & 
  theme(legend.position = 'bottom') 

############### size is problem 
ggsave(argv$fig6D, width = 22, height = 10, units = "cm") 

