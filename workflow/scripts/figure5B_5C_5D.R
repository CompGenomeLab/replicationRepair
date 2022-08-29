#### Packages and Libraries ####

library(argparser)
library(stringr)
library(ggpubr)
library(ggplot2)
library(dplyr)
library(reshape2)
library(patchwork)
set.seed(1) 

######## Arguments ##########
p <- arg_parser("producing the figure 5B,C,D")
p <- add_argument(p, "--noverlap", help="non-overlapping HeLa initiation zones")
p <- add_argument(p, "--overlap", help="overlapping HeLa initiation zones")
p <- add_argument(p, "--hela", help="mutation counts on overlapping HeLa initiation zones with scores")
p <- add_argument(p, "--noverlap_mut", help="mutation counts on non-overlapping HeLa initiation zones")
p <- add_argument(p, "--overlap_mut", help="mutation counts on overlapping HeLa initiation zones")
p <- add_argument(p, "--noverlap_TC", help="TC counts on non-overlapping HeLa initiation zones")
p <- add_argument(p, "--overlap_TC", help="TC counts on overlapping HeLa initiation zones")
p <- add_argument(p, "--data_prefix", help="name prefix of the dataframes that generate the plots")
p <- add_argument(p, "--fig5", help="figure output")

# Parse the command line arguments
argv <- parse_args(p)

hela2gm_imr <- argv$overlap

hela_no_overlap <- argv$noverlap

hela_bed <- argv$hela

hela2gmimr_bed <- argv$overlap_mut

helanover_bed <- argv$noverlap_mut
  
hela2gmimr_nuc <- argv$overlap_TC
  
helanover_nuc <- argv$noverlap_TC

#### Default Plot Format ####

source("workflow/scripts/plot_format.R")


#### Functions ####

source("workflow/scripts/functions.R")

#### plot A #####
# plot A will be a drawing
# we are creating an empty text for that part
p.A <- wrap_elements(grid::textGrob(''))

#### plot B ####

hela_overlap <- read.table( hela2gm_imr )
hela_no_overlap <- read.table( hela_no_overlap )
hela_overlap$V6 <- "overlapping"
hela_no_overlap$V6 <- "no overlap"

df <- rbind(hela_overlap, hela_no_overlap)
df$logV5 <- log2(df$V5)

write.table(df, file = paste0(argv$data_prefix, "B.csv"), quote = FALSE, 
            row.names = FALSE, sep = ",")

p.B <- ggboxplot(df, x = "V6", y = "logV5", outlier.shape = NA) + 
  stat_compare_means(comparisons = list(c(1,2)), label.y = 12.5)  +
  xlab("") +
  ylim(7.5, 13.5) +
  ylab("Initiation Zone\nScores (log2)") 

#### plot C ####

mut_df <- read.table( hela_bed )

mut_df <- mut_df[,c(1,2,3,4,5,9,11)]

#gm_df <- read.table( gm_bed )

mut_df$V11 <- as.numeric(mut_df$V11)

mut_df$windows <- as.numeric(
  as.character(str_split_fixed(mut_df$V4, "_",4)[,4])) - 101
mut_df$score <- as.numeric(
  as.character(str_split_fixed(mut_df$V4, "_",4)[,3]))
mut_df$domains <- as.character(str_split_fixed(mut_df$V4, "_",4)[,2])
mut_df$name <- as.character(str_split_fixed(mut_df$V4, "_",4)[,1])

mut_df_quant <- filter(mut_df, score != "Na", score != "NaN")
mut_df_quant <- within(mut_df_quant, score_quant <- as.integer(cut(score, quantile(score, probs=0:4/4), include.lowest=TRUE)))

mut_df_quant$domains <- factor(mut_df_quant$domains, levels = c("UTZ", "ERD", "DTZ", "LRD"))
mut_df_quant$V9 <- factor(
  mut_df_quant$V9, levels = c("+","-"))

mut_filt <- mut_df_quant[,c("name","domains","windows","score_quant","V9","V11")]

mut_agg <- aggregate(mut_filt$V11, 
                     by = list(mut_filt$name, 
                               mut_filt$domains, 
                               mut_filt$windows, 
                               mut_filt$score_quant,
                               mut_filt$V9),
                     FUN = mean)

mut_casted <- dcast(mut_agg, Group.1 + Group.2 + Group.3 + Group.4 ~ Group.5, 
                    value.var = "x")
mut_casted$plus_min <- mut_casted$"+" - mut_casted$"-" 

mut_plus <- filter(mut_casted, Group.3 > 0)
mut_plus_agg <- aggregate(x = mut_plus$plus_min, by = 
                            list(mut_plus$Group.1, mut_plus$Group.2, mut_plus$Group.4), 
                          FUN = "mean")
mut_plus_agg$direction <- "Right Replicating"
mut_minus <- filter(mut_casted, Group.3 < 0)
mut_minus_agg <- aggregate(x = mut_minus$plus_min, by = 
                             list(mut_minus$Group.1, mut_minus$Group.2, mut_minus$Group.4), 
                           FUN = "mean")
mut_minus_agg$direction <- "Left Replicating"
mut_data <- rbind(mut_plus_agg, mut_minus_agg) 

mut_data$Group.3 <- as.character(mut_data$Group.3)

write.table(mut_data, file = paste0(argv$data_prefix, "C.csv"), quote = FALSE, 
            row.names = FALSE, sep = ",")

p.C <- ggplot(data = mut_data, aes(fill = Group.3, x = Group.2, y = x)) + 
  geom_bar(position="dodge", stat = "identity") +
  facet_wrap(~direction) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  xlab("Replication Domains") + ylab(expression('raw MC'[p] - 'MC'[m]))

# adding and overriding the default plot format
p.C <- p.C + p_format + 
  labs(fill = "IZ Scores") +
  theme(legend.position = "right") +
  scale_fill_manual(values = c("gray80", "gray60", "gray40", "black"))

#### plot D ####

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
pD3_data <- rbind(mut_plus_agg, mut_minus_agg) 

write.table(mut, file = paste0(argv$data_prefix, "D1.csv"), quote = FALSE, 
            row.names = FALSE, sep = ",")

# create the plot
p.D.1 <- ggplot(mut) + 
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
p.D.1 <- p.D.1 + p_format + 
  theme(legend.position = "right",
        panel.border=element_rect(fill = NA),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.text.x = element_text(hjust=c(0.1, 0.5, 0.9))) 

mut$direction[mut$window_number<0] <- "Left Replicating" 
mut$direction[mut$window_number>0] <- "Right Replicating" 

mut <- filter(mut, direction != "NA")
mut$log2val <- log2(mut$norm)

pD1_boxplot <- mut[,c("repdomains", "log2val", "norm", "strands", "direction", "window_number")]

write.table(pD1_boxplot, file = paste0(argv$data_prefix, "D2.csv"), quote = FALSE, 
            row.names = FALSE, sep = ",")

p.D.2 <- ggplot(data=pD1_boxplot, aes(x=repdomains, y=norm, fill=strands)) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~direction) +
  xlab("Replication Domains") + 
  ylab("") +
  scale_fill_manual(values = c("+" = "#0571b0", 
                               "-" = "#ca0020"), 
                    guide = "none") +
  #scale_y_continuous(breaks = c(0, 1, 2),
  #                   limits = c(0, 2)) +
  labs(color = "Strands", fill = "") 

# adding and overriding the default plot format
p.D.2 <- p.D.2 + p_format + 
  stat_compare_means(label = "p.signif",  paired = TRUE, label.y = 0.037) + 
  theme(strip.background = element_blank(),
        legend.position = "right") 

#### Combining Plots with Patchwork ####

layout <- "
AABB
AACC
DDDD
DDDD
"

p.D <- p.D.1 + (p.D.2 + plot_layout(tag_level = 'new')) + 
  plot_layout(guides = "collect") + theme(legend.position = "bottom")

p.A + p.B + p.C + p.D +
  plot_annotation(tag_levels = 'A') +
  plot_layout(design = layout) & 
  theme(plot.tag = element_text(hjust = -0.2, size = 12, face="bold"))

ggsave(argv$fig5, width = 22, height = 18, units = "cm") 
