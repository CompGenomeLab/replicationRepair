#### Packages and Libraries ####

#library(DescTools)
library(dplyr)
library(ggplot2)
library(patchwork)
library(argparser)
set.seed(1) 

######## Arguments ##########
p <- arg_parser("producing the supplementary figure 5")
p <- add_argument(p, "--df", help="region file with read counts")
p <- add_argument(p, "--data_prefix", help="name prefix of the dataframes that generate the plots")
p <- add_argument(p, "-o", help="output")

# Parse the command line arguments
argv <- parse_args(p)

sample <- argv$df


#### Default Plot Format ####

source("workflow/scripts/plot_format.R")


#### Fuctions ####

source("workflow/scripts/functions.R")


#### Main ####

sample_df <- read.csv( sample, header = T )

sample_df$num <- 1
sample_df$eff[sample_df$ear_la > 0] <- 1
sample_df$eff[sample_df$ear_la < 0] <- 0
sample_df$eff <- as.numeric(sample_df$eff)

levels(sample_df$state_short) <- levels(as.factor(sample_df$state_short))

#tab <- xtabs(num ~ dataset_strand + eff + state_short, sample_df)
#Desc(tab)
#BreslowDayTest(tab)

#state <- filter(sample_df, state_short == "Transcription Associated")
#tab_state <- xtabs(num ~ dataset_strand + eff, state)
#Desc(tab_state)

df <- data.frame(yAxis = length(levels(sample_df$state_short)):1, 
                 boxOdds = log2(c(1.059, 1.071, 1.442, 1.143, 4.510, 2.6471, 
                                  2.215, 1.860, 1.529, 1.671)), 
                 boxCILow = log2(c(0.683, 0.595, 1.110, 0.401, 4.202, 0.2481, 
                                   2.160, 1.657, 1.003, 1.421)), 
                 boxCIHigh = log2(c(1.643, 1.925, 1.875, 3.258, 4.842, 28.2397, 
                                    2.272, 2.089, 2.332, 1.965))
)

df$yAxis <- levels(sample_df$state_short)
df$yAxis <- factor(
  df$yAxis, levels = rev(levels(sample_df$state_short)))

write.table(df, file = paste0(argv$data_prefix, ".csv"), 
            quote = FALSE, row.names = FALSE, sep = ",")

p <- ggplot(df, aes(x = boxOdds, y = yAxis)) + 
    geom_vline(aes(xintercept = 0), size = .25, linetype = "dashed") + 
    geom_errorbarh(aes(xmax = boxCIHigh, xmin = boxCILow), size = .5, height = 
                     .2, color = "gray50") +
    geom_point(size = 3.5, color = "orange") +
    theme_bw()+
    theme(panel.grid.minor = element_blank()) +
    ylab("") +
    xlab("log2 transformed Odds Ratio") +
    annotate(geom = "text", y =1.1, x = -5, 
           label = "Breslow-Day Test on Homogeneity of Odds Ratios\n                    (with Tarone correction)\n            p-value < 0.00000000000000022", 
             size = 3.5, hjust = 0) 

ggsave(argv$o, width = 22, height = 18, units = "cm") 

#sample_df$dataset_strand <- factor(
#  sample_df$dataset_strand, levels = c("LRD","ERD"))
#lap <- filter(sample_df, state_short == "Active Promoter")
#tab_lap <- xtabs(num ~ dataset_strand + eff, lap)
#Desc(tab_lap)
