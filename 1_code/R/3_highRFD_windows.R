#### library ####

library(ggplot2)
library(stringr)
library(dplyr)

#### prepare data ####

sourcePath <- paste("~/Documents/My_Projects/Project_Repair_Replication/", 
                    "Scripts/R/", sep = "") 

date <- "[2019.10.31]"

fr_name <- "HiRFD_windows"

source(paste(sourcePath, "2_report_sub_dfs.R", sep = ""))

#### set df, variables and rearrange ####

df <- fr_xr_ds_min_plus 

for (rearrange in 1) {  
  
  #### set variables ####
  
  dateout <- Sys.Date()
  dateout <- format(dateout, format = "[%Y.%m.%d]")

  window_number <- 11
  if (window_number %% 2 == 0) {
    half_window <- window_number / 2
  } else {
    half_window <- (window_number - 1) / 2 + 1
  }
  
  #### rearrange ####
  
  # window numbers separated from dataset names
  window <- data.frame(str_split_fixed(df$dataset, "_", -1))
  
  df$windows <- as.numeric(levels(
    window[ , ncol(window)]))[window[ , ncol(window)]] - half_window
  
  df$dataset <- "highRFD"
  
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
                    "Results/HiRFD/", sep = "")

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

figureInfo <- "3_highRFD_windows.R"

source(paste(sourcePath, "4_figure_info.R", sep = ""))



