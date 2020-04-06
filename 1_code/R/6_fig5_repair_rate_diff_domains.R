#### erd early late plot ####

library(ggplot2)
library(dplyr)

#### prepare data

sourcePath <- paste("~/Documents/myprojects/replicationRepair/1_code/R/", 
                    sep = "") 

date <- "[2019.10.31]"

fr_name <- "ERD_LRD_windows"

source(paste(sourcePath, "2_report_sub_dfs.R", sep = ""))

#### set df, variables and rearrange 

df = fr_xr_ds

for (rearrange in 1) {  
  
  #### set variables ####
  
  dateout <- Sys.Date()
  dateout <- format(dateout, format = "[%Y.%m.%d]")
  
  window_number <- 201
  if (window_number %% 2 == 0) {
    half_window <- window_number / 2
  } else {
    half_window <- (window_number - 1) / 2 + 1
  }
  
  #### rearrange ####
  
  # window numbers separated from dataset names
  
  df$windows <- as.numeric(gsub(".*_", "", df$dataset)) - half_window
  
  df$dataset <- gsub("_.*", "", df$dataset)
  
  rm(rearrange)
  
}

#### plot 

#### filter the data 

d <- filter(df, phase == "early" | phase == "late", dataset == "ERD", 
            replicate == "A")

#### add plot format

source(paste(sourcePath, "4_plot_format.R", sep = ""))

#### create the plot 

p.A.1.1 <- ggplot(d, aes(x = windows, y = log2(xr_ds))) + 
  geom_line(aes(linetype = sample_strand, color = phase)) +
  geom_hline(yintercept = 0, linetype="dashed", color="red") +
  facet_grid(~product~time_after_exposure~dataset,
             labeller = labeller(product = product_labs, 
                                 time_after_exposure = taex_labs)) +
  xlab(windows_lab) + ylab(fr_xr_ds_lab) +
  scale_x_continuous(limits = c(-half_window, half_window), 
                     breaks = c(-half_window, 0, half_window), 
                     labels = c("-1 Mbp", "0", "+1 Mbp")) + 
  scale_color_manual(name = "Phase", 
                     label = c("Early Phase", "Late Phase", "Asyncronized"), 
                     values = phase_colors) + 
  labs(color = "Phase", linetype = "Strands")

p.A.1.1 <- p.A.1.1 + p_format + ylim(-2, 2) +
  theme(strip.text.y = element_blank(),
        legend.position = "right")

p.A.1.1 # visualize

#### filter the data

d <- filter(df, phase == "early" | phase == "late", dataset == "LRD", 
            replicate == "A")

#### add plot format

source(paste(sourcePath, "4_plot_format.R", sep = ""))

#### create the plot

p.A.1.2 <- ggplot(d, aes(x = windows, y = log2(xr_ds))) + 
  geom_line(aes(linetype = sample_strand, color = phase)) +
  geom_hline(yintercept = 0, linetype="dashed", color="red") +
  facet_grid(~product~time_after_exposure~dataset,
             labeller = labeller(product = product_labs, 
                                 time_after_exposure = taex_labs)) +
  xlab(windows_lab) + ylab(fr_xr_ds_lab) +
  scale_x_continuous(limits = c(-half_window, half_window), 
                     breaks = c(-half_window, 0, half_window), 
                     labels = c("-1 Mbp", "0", "+1 Mbp")) + 
  scale_color_manual(name = "Phase", 
                     label = c("Early Phase", "Late Phase", "Asyncronized"), 
                     values = phase_colors) + 
  labs(color = "Phase", linetype = "Strands")

p.A.1.2 <- p.A.1.2 + p_format + ylim(-2, 2) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        strip.text.y = element_blank(),
        legend.position = "right")

p.A.1.2 # visualize

#### chromosomes diff ####

library(ggplot2)
library(stringr)
library(ggpubr)

#### prepare data

sourcePath <- paste("~/Documents/myprojects/replicationRepair/1_code/R/", 
                    sep = "") 

date <- "[2019.10.31]"

fr_name <- "repdomain"

source(paste(sourcePath, "2_report_sub_dfs.R", sep = ""))

#### set df, variables and rearrange 

dateout <- Sys.Date()
dateout <- format(dateout, format = "[%Y.%m.%d]")
options(scipen = 999)

df <- fr_xr_ds

df$chromosomes <- factor(df$chromosomes, levels = c("chr1", "chr2", "chr3", 
                                                    "chr4", "chr5", "chr6", 
                                                    "chr7", "chr8", "chr9", 
                                                    "chr10", "chr11", "chr12", 
                                                    "chr13", "chr14", "chr15", 
                                                    "chr16", "chr17", "chr18", 
                                                    "chr19", "chr20", "chr21", 
                                                    "chr22", "chrX"))

#### plot 

#### filter the data 

d = filter(df, dataset == "ERD" | dataset == "LRD", replicate == "A", 
           product == "CPD", time_after_exposure == "12")

#### add plot format 

source(paste(sourcePath, "4_plot_format.R", sep = ""))

#### create the plot 

p.B.1.1 <- ggplot(d, aes(x = chromosomes, y = log2(xr_ds))) + 
  geom_boxplot(aes(fill = dataset), outlier.shape = NA) +
  geom_hline(yintercept = -1, linetype="dashed", color = "darkred") +
  geom_hline(yintercept = 1, linetype="dashed", color = "darkgreen") +
  facet_grid(phase~time_after_exposure~product,
             labeller = labeller(product = product_labs, 
                                 time_after_exposure = taex_labs, 
                                 replicate = rep_labs, phase = phase_labs)) +
  xlab("Chromosomes") + ylab(fr_xr_ds_lab) +
  scale_fill_manual(name="", 
                    values = repdomain_colors) +
  ylim(-3,2)


p.B.1.1 <- p.B.1.1 +  p_format + theme(axis.text.x = 
                             element_text(size = 12, vjust = 0.6, angle = 60))

p.B.1.1 # visualize

#### samples diff ####

library(ggplot2)
library(stringr)
library(ggpubr)

#### prepare data 

sourcePath <- paste("~/Documents/myprojects/replicationRepair/1_code/R/", 
                    sep = "") 

date <- "[2019.10.31]"

fr_name <- "repdomain"

source(paste(sourcePath, "2_report_sub_dfs.R", sep = ""))

#### set df, variables and rearrange 

dateout <- Sys.Date()
dateout <- format(dateout, format = "[%Y.%m.%d]")
options(scipen = 999)

df <- fr_xr_ds

#### plot

#### filter the data 

d = filter(df, dataset == "ERD" | dataset == "LRD", replicate == "A", 
           phase != "async")

#### add plot format 

source(paste(sourcePath, "4_plot_format.R", sep = ""))

phase_labs <- c("Asyncronized", "Early Phase", "Late Phase")
names(phase_labs) <- c("async", "early", "late")

#### create the plot 

p.A.1.3 <- ggplot(d, aes(x = dataset, y = log2(xr_ds))) + 
  geom_boxplot(aes(fill = dataset), outlier.shape = NA) +
  geom_hline(yintercept = 0, linetype="dashed", color="red") +
  facet_grid(product~time_after_exposure~phase, 
             labeller = labeller(product = product_labs, 
                                 time_after_exposure = taex_labs, 
                                 phase = phase_labs, replicate = rep_labs)) + 
  xlab(repdo_lab) + ylab(fr_xr_ds_lab) +
  scale_fill_manual(name="", 
                    values = repdomain_colors) +
  guides(fill = FALSE) +
  ylim(-2,2) 


p.A.1.3 <- p.A.1.3 +  p_format +  stat_compare_means(method = "t.test", 
                                         label.x = 1.3, label.y = 1.8) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

p.A.1.3 # visualize

#### chr states ####

library(ggplot2)
library(stringr)
library(dplyr)

#### prepare data

sourcePath <- paste("~/Documents/myprojects/replicationRepair/1_code/R/", 
                    sep = "") 

date <- "[2019.11.07]"

fr_name <- "ChrStates_chromhmm_windows"

source(paste(sourcePath, "2_report_sub_dfs.R", sep = ""))

#### set df, variables and rearrange

dateout <- Sys.Date()
dateout <- format(dateout, format = "[%Y.%m.%d]")

df <- fr_xr_ds_ear_la

for (rearrange in 1) {  
  
  #### set variables 
  
  chr_states <- c("Tss", "TssF", "PromF", "PromP", "Enh", "EnhF", "EnhWF", 
                  "EnhW", "FaireW", "Ctcf", "Gen5'", "Elon", "ElonW", "Gen3'", 
                  "Pol2", "H4K20", "Low", "ReprD", "Repr", "ReprW", 
                  "Quies", "Art")
  
  state <- c("Active Promoter", "Active Promoter", "Promoter Flanking", 
             "Inactive Promoter", "Candidate Strong Enhancer", 
             "Candidate Strong Enhancer", "Candidate Weak Enhancer", 
             "Candidate Weak Enhancer", "Candidate Weak Enhancer", 
             "Distal CTCF/Candidate Insulator", "Transcription Associated", 
             "Transcription Associated", "Transcription Associated", 
             "Transcription Associated", "Transcription Associated", 
             "Transcription Associated", 
             "Low Activity Proximal to Active States", "Polycob Repressed", 
             "Polycob Repressed", "Polycob Repressed", 
             "Heterochromatin/Repetitive/ \n  Copy Number Variation", 
             "Heterochromatin/Repetitive/ \n  Copy Number Variation")
  
  #### rearrange 
  
  df$dataset <- factor(df$dataset, levels = c("Tss", "TssF", "PromF", "PromP", 
                                              "Enh", "EnhF", "EnhWF", "EnhW", 
                                              "FaireW", "Ctcf", "Gen5'", 
                                              "Elon", "ElonW", "Gen3'", "Pol2", 
                                              "H4K20", "Low", "ReprD", "Repr", 
                                              "ReprW", "Quies", "Art"))
  
  df$states <- NA
  
  for (i in 1:length(chr_states)) {
    
    df <- within(df, states[dataset == chr_states[i]] <- state[i])
    
  }
  
  df$states <- as.factor(df$states)
  
  df$states <- factor(
    df$states, levels = c("Active Promoter", 
                          "Promoter Flanking", 
                          "Inactive Promoter", 
                          "Candidate Strong Enhancer", 
                          "Candidate Weak Enhancer", 
                          "Distal CTCF/Candidate Insulator", 
                          "Transcription Associated", 
                          "Low Activity Proximal to Active States", 
                          "Polycob Repressed", 
                          "Heterochromatin/Repetitive/ \n  Copy Number Variation"))
  
  rm(state, chr_states)
  
  
  
}

#### plot

#### filter the data

df = filter(df, ear_la != "NaN" & ear_la != "Inf" )

d = filter(df, dataset_strand != "DTZ" & dataset_strand != "UTZ" & 
             product == "CPD" & replicate == "A", 
           time_after_exposure == "12")

d = aggregate(x = d$ear_la, by = 
                list(d$dataset, d$dataset_strand, d$states), FUN = "mean")

#### add plot format  

source(paste(sourcePath, "4_plot_format.R", sep = ""))

#### create the plot 

p.B.1.2 <- ggplot(d, aes(x = Group.1, y = log2(x))) + 
  geom_hline(yintercept = 0.2, linetype = "dashed", color = "#7fcdbb", 
             size = 2) +
  geom_hline(yintercept = -0.2, linetype = "dashed", color = "#2c7fb8",
             size = 2) +
  geom_point(shape = 21, size = 4, aes(fill = factor(Group.3))) +
  facet_wrap(~Group.2) +
  xlab(chrState_lab) + ylab(fr_xr_ds_ear_la_lab) +
  scale_fill_manual(name = "Chromatin States", values = state_colors) +
  guides(size = FALSE) 

p.B.1.2 <- p.B.1.2 + p_format + theme(legend.position = "right", 
                          axis.text.x = element_text(size = 12, vjust = 0.6, 
                                                     angle = 60))

p.B.1.2 # visualize


#### final plot ####

((p.A.1.1 + p.A.1.2 + p.A.1.3) + plot_layout(guides = "collect")) / 
  p.B.1.1 / p.B.1.2 + plot_layout(tag_level = 'new') +  
  plot_annotation(tag_levels = c('A'), tag_suffix = ':')

ggsave("~/Desktop/fig5.png", width = 497, height = 410, units = "mm")




