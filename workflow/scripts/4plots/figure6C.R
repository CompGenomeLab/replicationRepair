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
hela_bed <- "/Users/azgarian/Documents/myprojects/replicationRepair/results/processed_files/processed_files/[2021.04.05.10.47]iz_hela_to_gm_imr_repdomains_with_scores_201_100_melanoma_mutations_combined_intergenic.bed"


#### Default Plot Format ####

source("/Users/azgarian/Documents/myprojects/replicationRepair/workflow/scripts/plot_format.R")


#### Functions ####

source("/Users/azgarian/Documents/myprojects/replicationRepair/workflow/scripts/functions.R")


#### main ####

mut_df <- read.table( hela_bed )

#gm_df <- read.table( gm_bed )

mut_df$V7 <- as.numeric(mut_df$V7)

mut_df$windows <- as.numeric(
  as.character(str_split_fixed(mut_df$V4, "_",4)[,4])) - 101
mut_df$score <- as.numeric(
  as.character(str_split_fixed(mut_df$V4, "_",4)[,3]))
mut_df$domains <- as.character(str_split_fixed(mut_df$V4, "_",4)[,2])
mut_df$name <- as.character(str_split_fixed(mut_df$V4, "_",4)[,1])

mut_df_quant <- filter(mut_df, score != "Na", score != "NaN")
mut_df_quant <- within(mut_df_quant, score_quant <- as.integer(cut(score, quantile(score, probs=0:4/4), include.lowest=TRUE)))



mut_df_quant$domains <- factor(mut_df_quant$domains, levels = c("UTZ", "ERD", "DTZ", "LRD"))
mut_df_quant$V6 <- factor(
  mut_df_quant$V6, levels = c("+","-"))


# for plot A.3

mut_filt <- mut_df_quant[,c("name","domains","windows","score_quant","V6","V7")]

mut_agg <- aggregate(mut_filt$V7, 
                     by = list(mut_filt$name, 
                               mut_filt$domains, 
                               mut_filt$windows, 
                               mut_filt$score_quant,
                               mut_filt$V6),
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

p.A <- ggplot(data = mut_data, aes(fill = Group.3, x = Group.2, y = x)) + 
  geom_bar(position="dodge", stat = "identity") +
  facet_wrap(~direction) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  xlab("Replication Domains") + ylab(expression('raw MC'[p] - 'MC'[m]))

# adding and overriding the default plot format
p.A <- p.A + p_format + 
  labs(fill = "IZ Scores") +
  theme(legend.position = "right") +
  scale_fill_manual(values = c("gray80", "gray60", "gray40", "black"))

p.A


ggsave("~/Desktop/fig6C.png", width = 12, height = 4.5, units = "cm")
