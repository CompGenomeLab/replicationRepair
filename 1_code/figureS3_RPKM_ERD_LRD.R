#### Packages and Libraries ####

library(stringr)
library(ggpubr)
library(ggplot2)
library(dplyr)
library(reshape2)
library(patchwork)
library(ggthemes)


#### Variables ####

# name of the sample csv file for plot A
sample_csv_pA <- paste("~/Documents/myprojects/replicationRepair/3_output/",
                       "gitignore/1_TextforPlotting/", 
                       "[2019.10.31]final_report_ERD_LRD_windows_ready.csv", 
                       sep = "")

# path of the default plot format and functions
sourcePath <- "~/Documents/myprojects/replicationRepair/1_code/r/"


#### Default Plot Format ####

source(paste(sourcePath, "4_plot_format.R", sep = ""))


#### Fuctions ####

source(paste(sourcePath, "4_functions.R", sep = ""))


#### Main ####

# for plot A
pA_sample_df <- read.csv( sample_csv_pA )
pA_df_org <- window_numbering( pA_sample_df, 4, 101 )
pA_df_org$dataset <- gsub("_.*", "", pA_df_org$dataset)

# filtering samples 
pA1_data <- filter(pA_df_org, phase != "async", 
                   dataset == "ERD" | dataset == "LRD",
                   replicate == rep, time_after_exposure == "12",
                   product == "CPD")


#### Plot A.1 ####

# create the plot
p.A.1 <- ggplot(pA1_data, aes(x = windows, y = RPKM)) + 
  geom_smooth(aes(color = phase, linetype = sample_strand), se=FALSE) +
  facet_grid(~product~time_after_exposure~method~dataset,
             labeller = labeller(product = product_labs, 
                                 time_after_exposure = taex_labs,
                                 method = method_labs)) +
  xlab("Relative Position (kb)") + ylab(fr_lab) +
  scale_x_continuous(limits = c(-101, 101), 
                     breaks = c(-101, 0, 101), 
                     labels = c("-1000", "0", "+1000")) + 
  scale_color_manual(name = "Phase", 
                     label = c("Early Phase", "Late Phase", "Asyncronized"), 
                     values = phase_colors) + 
  labs(color = "Phase", linetype = "Strands") 

# adding and overriding the default plot format
p.A.1 <- p.A.1 + p_format 


#### Combining Plots with Patchwork ####


ggsave("~/Desktop/supfig3.svg", 
       width = 22, height = 18, units = "cm")

