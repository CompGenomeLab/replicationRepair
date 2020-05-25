#### chr states ####

library(ggpubr)
library(ggplot2)
library(stringr)
library(dplyr)
library(ggthemes)
library(patchwork)

#### prepare data

sourcePath <- paste("~/Documents/myprojects/replicationRepair/1_code/R/", 
                    sep = "") 

date <- "[2020.04.07]"

fr_name <- "chromhmm_windows_chr"

source(paste(sourcePath, "2_report_sub_dfs.R", sep = ""))

source(paste(sourcePath, "4_plot_format.R", sep = ""))

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
             product == "64_PP" & replicate == "B", 
           time_after_exposure == "12")

d = aggregate(x = d$ear_la, by = 
                list(d$dataset, d$dataset_strand, d$states, d$chromosomes), FUN = "mean")

#### add plot format  

p_format <- theme_pubr() +
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12, vjust = 0.6),
        axis.text.y = element_text(size = 12, vjust = 0.1),
        strip.text.x = element_text(size = 16),
        strip.text.y = element_text(size = 16), #angle = 360),
        strip.background = element_blank(),
        legend.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 16),
        legend.position = "bottom")

#### create the plot 

p.B.1.2 <- ggplot(d, aes(x = Group.1, y = log2(x))) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "blue") +
  geom_boxplot(aes(fill = factor(Group.3)), outlier.shape = NA) +
  #geom_jitter(aes(fill = factor(Group.3))) +
  facet_wrap(~Group.2) +
  xlab(chrState_lab) + ylab(fr_xr_ds_ear_la_lab) +
  scale_fill_manual(name = "Chromatin States", values = state_colors) +
  guides(size = FALSE) 

p.B.1.2 <- p.B.1.2 + p_format + ylim(-1,1) +
  #theme_bw() +
  theme(
    panel.spacing = unit(1.5, "lines"),
    panel.border=element_rect(size=1, fill = NA),
    strip.text = element_blank(),
    legend.position = "right", 
    axis.text.x = element_text(size = 12, vjust = 0.6, 
                               angle = 60))

#p.B.1.2 # visualize

#ggsave("~/Desktop/chrom2_64_12_repA.png", width = 510, height = 210, units = "mm")


#### repair rate ####

df <- fr_xr_ds

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

df = filter(df, xr_ds != "NaN" & xr_ds != "Inf" )

d = filter(df, dataset_strand != "DTZ" & dataset_strand != "UTZ" & 
             product == "64_PP" & replicate == "B" & phase != "async" & 
           time_after_exposure == "12")

d = aggregate(x = d$xr_ds, by = 
                list(d$dataset, d$dataset_strand, d$phase, d$states, d$chromosomes), FUN = "mean")

#### add plot format  

p_format <- theme_pubr() +
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12, vjust = 0.6),
        axis.text.y = element_text(size = 12, vjust = 0.1),
        strip.text.x = element_text(size = 16),
        strip.text.y = element_text(size = 16), #angle = 360),
        strip.background = element_blank(),
        legend.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 16),
        legend.position = "bottom")

#### create the plot 

p.B.1.1 <- ggplot(d, aes(x = Group.1, y = log2(x))) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_boxplot(aes(fill = factor(Group.4)), outlier.shape = NA) +
  #geom_jitter(aes(fill = factor(Group.4))) +
  facet_grid(~Group.3~Group.2, labeller = labeller(Group.3 = phase_labs)) +
  xlab("") + ylab(fr_xr_ds_lab) +
  scale_fill_manual(name = "Chromatin States", values = state_colors) +
  guides(size = FALSE) 

p.B.1.1 <- p.B.1.1 + p_format + ylim(-2,3) +
  #theme_bw() +
  theme(
    panel.spacing = unit(1.5, "lines"),
    panel.border=element_rect(size=1, fill = NA),
    #strip.text = element_blank(),
    legend.position = "right", 
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank())

#p.B.1.1 # visualize

# final plot

p.B.1.1 / p.B.1.2 + plot_layout(guides = "collect") 

ggsave("~/Desktop/rr_chrom2_64PP_12_repB.png", width = 510, height = 310, units = "mm")
