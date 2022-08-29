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
                 boxOdds = log2(c(1.521, 1.880, 0.3684, 0.826, 1.427, 1.046,
                                 1.713, 2.251, 1.677, 4.561)), 
                 boxCILow = log2(c(0.973, 1.272, 0.0353, 0.510, 1.123, 0.443, 
                              1.474, 2.195, 1.499, 4.251)), 
                 boxCIHigh = log2(c(2.375, 2.779, 3.8504, 1.340, 1.813, 2.469, 
                               1.991, 2.308, 1.875, 4.894))
)

df$yAxis <- levels(sample_df$state_short)
df$yAxis <- factor(
  df$yAxis, levels = rev(levels(sample_df$state_short)))

df_erd <- data.frame(yAxis = length(levels(sample_df$state_short)):1, 
                 boxOdds = log2(c(1.188, 1.358, 0.5263, 0.908, 1.186, 1.020,
                             1.288, 1.483, 1.285, 2.055)), 
                 boxCILow = log2(c(1.010, 1.149, 0.0935, 0.703, 1.067, 0.701, 
                              1.211, 1.465, 1.211, 1.969)), 
                 boxCIHigh = log2(c(1.397, 1.605, 2.9639, 1.172, 1.319, 1.484, 
                               1.369, 1.501, 1.363, 2.144))
)

df_erd$yAxis <- levels(sample_df$state_short)
df_erd$yAxis <- factor(
  df_erd$yAxis, levels = rev(levels(sample_df$state_short)))

df_lrd <- data.frame(yAxis = length(levels(sample_df$state_short)):1, 
                     boxOdds = log2(c(0.781, 0.722, 1.4286, 1.099, 0.831, 0.975,
                                 0.752, 0.659, 0.766, 0.450)), 
                     boxCILow = log2(c(0.588, 0.576, 0.7548, 0.875, 0.728, 0.601, 
                                  0.688, 0.650, 0.727, 0.437)), 
                     boxCIHigh = log2(c(1.038, 0.905, 2.7038, 1.381, 0.950, 1.582, 
                                   0.822, 0.668, 0.808, 0.465))
)

df_lrd$yAxis <- levels(sample_df$state_short)
df_lrd$yAxis <- factor(
  df_lrd$yAxis, levels = rev(levels(sample_df$state_short)))

write.table(df, file = paste0(argv$data_prefix, ".csv"), 
            quote = FALSE, row.names = FALSE, sep = ",")
#write.table(df_erd, file = paste0(argv$data_prefix, "B1.csv"), 
#            quote = FALSE, row.names = FALSE, sep = ",")
#write.table(df_lrd, file = paste0(argv$data_prefix, "B2.csv"), 
#            quote = FALSE, row.names = FALSE, sep = ",")

p_all <- ggplot(df, aes(x = boxOdds, y = yAxis)) + 
    geom_vline(aes(xintercept = 0), size = .25, linetype = "dashed") + 
    geom_errorbarh(aes(xmax = boxCIHigh, xmin = boxCILow), size = .5, height = 
                     .2, color = "gray50") +
    geom_point(size = 3.5, color = "orange") +
    theme_bw()+
    theme(panel.grid.minor = element_blank()) +
    ylab("") +
    xlab("log2 transformed Odds Ratio") +
    annotate(geom = "text", y =1.1, x = -5, 
             label = "Breslow-Day test on Homogeneity of Odds Ratios\n          p-value < 0.00000000000000022", 
             size = 3.5, hjust = 0) 

p_rel <- ggplot(df, aes(x = boxOdds, y = yAxis)) + 
  geom_vline(aes(xintercept = 0), size = .25, linetype = "dashed") + 
  geom_errorbarh(data=df_erd, aes(y = yAxis, 
                   xmax = boxCIHigh, xmin = boxCILow), size = .5, height = 
                   .2, color = "red") +
  geom_point(data=df_erd, aes(x = boxOdds, y = yAxis), 
             size = 3.5, color = "red") +
  geom_errorbarh(data=df_lrd, aes(y = yAxis, 
                                  xmax = boxCIHigh, xmin = boxCILow), size = .5, height = 
                   .2, color = "blue") +
  geom_point(data=df_lrd, aes(x = boxOdds, y = yAxis), 
             size = 3.5, color = "blue") +
  theme_bw()+
  theme(panel.grid.minor = element_blank()) +
  ylab("") +
  xlab("log2 transformed Relative Risks") 

p_all

ggsave(argv$o, width = 22, height = 18, units = "cm") 

#sample_df$dataset_strand <- factor(
#  sample_df$dataset_strand, levels = c("LRD","ERD"))
#lap <- filter(sample_df, state_short == "Active Promoter")
#tab_lap <- xtabs(num ~ dataset_strand + eff, lap)
#Desc(tab_lap)
