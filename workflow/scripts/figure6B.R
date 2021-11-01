#### Packages and Libraries ####

library(argparser)
library(stringr)
library(ggpubr)
library(ggplot2)
library(dplyr)
library(reshape2)
library(patchwork)

######## Arguments ##########
p <- arg_parser("producing the figure 6B")
p <- add_argument(p, "--noverlap", help="non-overlapping HeLa initiation zones")
p <- add_argument(p, "--overlap", help="overlapping HeLa initiation zones")
p <- add_argument(p, "--fig6B", help="figure output")

# Parse the command line arguments
argv <- parse_args(p)

hela2gm_imr <- argv$overlap

hela_no_overlap <- argv$noverlap

#### Default Plot Format ####

source("workflow/scripts/plot_format.R")


#### Functions ####

source("workflow/scripts/functions.R")


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

ggsave(argv$fig6B, width = 8, height = 4.5, units = "cm")
