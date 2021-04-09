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
hela_bed <- "/home/azgarian/Desktop/repairRep_revision/inZones/okseq/hela2gm.txt"
  
gm_bed <- "/home/azgarian/Desktop/repairRep_revision/inZones/okseq/gm2hela.txt"

hela2gm_imr <- "/home/azgarian/Desktop/repairRep_revision/inZones/hela_intersect2_gm_imr.bed"

hela_no_overlap <- "/home/azgarian/Desktop/repairRep_revision/inZones/hela_no_overlap.bed"

# path of the default plot format and functions
sourcePath <- "~/Documents/myprojects/replicationRepair/1_code/r/"


#### Default Plot Format ####

source(paste(sourcePath, "4_plot_format.R", sep = ""))


#### Fuctions ####

source(paste(sourcePath, "4_functions.R", sep = ""))


#### main ####

hela_df <- read.table( hela_bed )

gm_df <- read.table( gm_bed )

hela_df <- within(hela_df, V6[V7 == '-1'] <- 'no overlap')
hela_df <- within(hela_df, V6[V7 != '-1'] <- 'overlapping')
hela_df$name <- "hela"

gm_df <- within(gm_df, V6[V7 == '-1'] <- 'no overlap')
gm_df <- within(gm_df, V6[V7 != '-1'] <- 'overlapping')
gm_df$name <- "gm06990"

df <- rbind(hela_df, gm_df)

df_agg <- df[,c("V4", "V6", "name")]
df_agg <- aggregate(V4 ~ V6 + name, data = df_agg, mean)

p <- ggplot(data = df_agg, aes(x = V6, y = V4)) + 
  geom_bar(stat = "identity") +
  facet_grid(~name) +
  xlab("") + ylab("Initiation Zone Scores") 

# adding and overriding the default plot format
p <- p + p_format 

p

ggsave("~/Desktop/hela_gm_intersect_barplot.png", width = 22, height = 9, units = "cm") 

df$logV4 <- log2(df$V4)
df$name2 <- paste(df$V6, df$name)
df$name2 <- factor(df$name2, levels = c("overlapping hela", "no overlap hela", "overlapping gm06990", "no overlap gm06990"))

ggboxplot(df, x = "name2", y = "logV4") + 
  stat_compare_means(comparisons = list(c(1,2), c(3,4)), label.y = 16)  +
  xlab("") +
  ylab("Initiation Zone Scores (log2)")

ggsave("~/Desktop/hela_gm_intersect_boxplot.png", width = 22, height = 18, units = "cm")






hela_overlap <- read.table( hela2gm_imr )
hela_no_overlap <- read.table( hela_no_overlap )
hela_overlap$V6 <- "overlapping"
hela_no_overlap$V6 <- "no overlap"

df <- rbind(hela_overlap, hela_no_overlap)
df$logV5 <- log2(df$V5)

ggboxplot(df, x = "V6", y = "logV5", outlier.shape = NA) + 
  stat_compare_means(comparisons = list(c(1,2)), label.y = 12.5)  +
  xlab("") +
  ylim(7.5, 13.5) +
  ylab("Initiation Zone\nScores (log2)") 

ggsave("~/Desktop/hela_intersect_gm_imr90_boxplot.png", width = 8, height = 4.5, units = "cm")
