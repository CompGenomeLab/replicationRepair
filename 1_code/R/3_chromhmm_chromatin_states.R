#### library ####

library(ggplot2)
library(stringr)
library(dplyr)

#### prepare data ####

sourcePath <- paste("~/Documents/My_Projects/Project_Repair_Replication/", 
                    "Scripts/R/", sep = "") 

date <- "[2019.11.07]"

fr_name <- "ChrStates_chromhmm_windows"

source(paste(sourcePath, "2_report_sub_dfs.R", sep = ""))

#### set df, variables and rearrange ####

dateout <- Sys.Date()
dateout <- format(dateout, format = "[%Y.%m.%d]")

setwd(paste("~/Documents/My_Projects/Project_Repair_Replication/", 
            "Results/Chromatin_states/", sep = ""))

df <- fr_xr_ds_ear_la

for (rearrange in 1) {  
  
  #### set variables ####
  
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
  
  #### rearrange ####
  
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

#### plot ####

#### filter the data ####

df = filter(df, ear_la != "NaN" & ear_la != "Inf" )

d = filter(df, dataset_strand != "DTZ" & dataset_strand != "UTZ" & 
             product == "CPD" & replicate == "A", 
           time_after_exposure == "12", sample_strand == "+")

#### add plot format #### 

source(paste(sourcePath, "4_plot_format.R", sep = ""))

#### create the plot ####

p <- ggplot(d, aes(x = dataset, y = log2(ear_la), color = dataset_strand)) + 
  geom_point(shape = 21, size = 4, stroke = 2, aes(fill = factor(states))) +
  geom_hline(yintercept = 0.2, linetype = "dashed", color = "black") +
  geom_hline(yintercept = -0.2, linetype = "dashed", color = "black") +
  xlab(chrState_lab) + ylab(fr_xr_ds_ear_la_lab) +
  annotate(geom = "label", x = 3.1, y = 0.21, 
           label = "Repaired Dominantly at Early Phase", 
           fill = "black", color = "white") +
  annotate(geom = "label", x = 3, y = -0.21, 
           label = "Repaired Dominantly at Late Phase", 
           fill = "black", color = "white") +
  scale_fill_manual(name = "Chromatin States", values = state_colors) +
  scale_color_manual(name = "Replication Domains", values = repdomain_colors) +
  guides(size = FALSE) 

p <- p + p_format + theme(legend.position = "right", 
                          axis.text.x = element_text(size = 12, vjust = 0.6, 
                                                     angle = 60))

p # visualize

#### save the plot ####

figurePath <- paste("~/Documents/My_Projects/Project_Repair_Replication/", 
                    "Results/Chromatin_states/", sep = "")

figureName <- paste(dateout, "Phase_Differences_of_Domains_Between_", 
                    "Chromatin_States.pdf", sep = "")

ggsave(path = figurePath, filename = figureName, 
       width = 397, height = 210, units = "mm")

figurePNG <- sub(".pdf", ".png", figureName)

ggsave(path = figurePath, filename = figurePNG, 
       width = 397, height = 210, units = "mm")

#### save the figure info ####

dataInfo <- paste(date, "final_report_", fr_name, "_info.TXT", sep = "")

figureInfo <- "3_chromhmm_chromatin_states.R"

source(paste(sourcePath, "4_figure_info.R", sep = ""))





