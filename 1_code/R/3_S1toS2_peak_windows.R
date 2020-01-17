#### library ####

library(ggplot2)
library(stringr)
library(dplyr)

#### prepare data ####

sourcePath <- paste("~/Documents/My_Projects/Project_Repair_Replication/", 
                    "Scripts/R/", sep = "") 

date <- "[2019.11.06]"

fr_name <- "S1toS2_windows"

source(paste(sourcePath, "2_report_sub_dfs.R", sep = ""))

#### set df, variables and rearrange ####

df <- fr_xr_ds

for (rearrange in 1) {  
 
  #### set variables ####
  
  dateout <- Sys.Date()
  dateout <- format(dateout, format = "[%Y.%m.%d]")

  window_number <- 101
  if (window_number %% 2 == 0) {
    half_window = window_number / 2
  } else {
    half_window = (window_number - 1) / 2 + 1
  }
  
  
  #### rearrange ####
  
  # new column direction added and filled with values "right_replicating" 
  # and "left_replicating"
  
  df$direction <- NA
  
  df <- within(df, direction[dataset_strand == "+"] <- "right_replicating")
  
  df <- within(df, direction[dataset_strand == "-"] <- "left_replicating")
  
  df$direction <- as.factor(df$direction)
  
  # new column strand added and filled with values "leading" and "lagging"
  df$strand <- NA
  
  df <- within(df, strand[(dataset_strand == "+" & sample_strand == "+") | 
                            (dataset_strand == "-" & sample_strand == "-")] 
               <- "lagging")
  
  df <- within(df, strand[(dataset_strand == "+" & sample_strand == "-") |
                            (dataset_strand == "-" & sample_strand == "+")] 
               <- "leading")
  
  df$strand <- as.factor(df$strand)
  
  # window numbers separated from dataset names
  window <- data.frame(str_split_fixed(df$dataset, "_", -1))
  
  df$windows <- as.numeric(levels(
    window[ , ncol(window)]))[window[ , ncol(window)]] - half_window
  
  df$dataset <- "S1toS2"
  
  rm(window, rearrange, half_window)

}

#### plot ####

#### filter the data ####

d <- filter(df, phase == "early" | phase == "late") # filter

#### add plot format #### 

source(paste(sourcePath, "4_plot_format.R", sep = ""))

#### create the plot ####

p <- ggplot(d) # plot

p <- p + p_format # add plot format

p # visualize

#### save the plot ####

figurePath <- paste("~/Documents/My_Projects/Project_Repair_Replication/", 
                    "Results/S1toS2_regions/", sep = "")

figureName <- paste(dateout, "", 
                    "", fr_name, 
                    "", sep = "")

ggsave(path = figurePath, filename = figureName, 
       width = 297, height = 210, units = "mm")

figurePNG <- sub(".pdf", ".png", figureName)

ggsave(path = figurePath, filename = figurePNG, 
       width = 297, height = 210, units = "mm")

#### save the figure info ####

dataInfo <- paste(date, "final_report_", fr_name, "_info.TXT", sep = "")

figureInfo <- "3_S1toS2_peak_windows.R"

source(paste(sourcePath, "4_figure_info.R", sep = ""))


