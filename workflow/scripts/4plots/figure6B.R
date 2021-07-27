#### Packages and Libraries ####

library(stringr)
library(ggpubr)
library(ggplot2)
library(dplyr)
library(reshape2)
library(patchwork)
library(ggthemes)


#### Variables ####
hela2gm_imr <- "/home/azgarian/Desktop/repairRep_revision/inZones/hela_intersect2_gm_imr.bed"

hela_no_overlap <- "/home/azgarian/Desktop/repairRep_revision/inZones/hela_no_overlap.bed"


#### Default Plot Format ####

source("4_plot_format.R")


#### Fuctions ####

source("4_functions.R")


#### main ####

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

ggsave("~/Desktop/fig6B.png", width = 8, height = 4.5, units = "cm")
