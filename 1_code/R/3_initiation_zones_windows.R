#### library ####

library(ggplot2)
library(stringr)
library(dplyr)

#### prepare data ####

sourcePath <- paste("~/Documents/My_Projects/Project_Repair_Replication/", 
                    "Scripts/R/", sep = "") 

date <- "[2019.10.31]"

fr_name <- "InZones_windows"

source(paste(sourcePath, "2_report_sub_dfs.R", sep = ""))

#### set df, variables and rearrange ####

df <- fr_xr_ds_min_plus

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
  window <- data.frame(str_split_fixed(df$dataset, "_", -1))
  
  df$windows <- as.numeric(levels(
    window[ , ncol(window)]))[window[ , ncol(window)]] - half_window
  
  df$dataset <- "Initiation Zones"
  
  rm(window, rearrange)
  
  
}

#### plot ####

#### filter the data ####

d <- filter(df, phase == "early" | phase == "late")

#### add plot format #### 

source(paste(sourcePath, "4_plot_format.R", sep = ""))

#### create the plot ####

p <- ggplot(d, aes(x = windows, y = log2(min_plus))) + 
  geom_line(aes(color = phase)) + 
  geom_line(y=0, color="red", linetype="dashed") +
  facet_grid(~product~time_after_exposure~replicate, 
             labeller = labeller(product = product_labs, 
                                 time_after_exposure = taex_labs, 
                                 replicate = rep_labs)) + 
  xlab(windows_lab) + ylab(fr_xr_ds_min_plus_lab) +
  scale_x_continuous(limits = c(-half_window-5, half_window+5), 
                     breaks = c(-half_window, 0, half_window), 
                     labels = c("-10 kb", "Initiation Zones", "+10 kbp")) + 
  scale_color_manual(name = "Phase", 
                     label = c("Early Phase", "Late Phase", "Asyncronized"), 
                     values = phase_colors) + 
  labs(color = "Strands", linetype = "Strands")

p <- p + p_format

p # visualize

#### save the plot ####

figurePath <- paste("~/Documents/My_Projects/Project_Repair_Replication/", 
                    "Results/InZones/", sep = "")

figureName <- paste(dateout, "Minus_over_Plus_Strand_Value_of_", 
                    "Repair_over_Damage_Value_of_", fr_name, 
                    "_for_Every_Sample.pdf", sep = "")

ggsave(path = figurePath, filename = figureName, 
       width = 297, height = 210, units = "mm")

figurePNG <- sub(".pdf", ".png", figureName)

ggsave(path = figurePath, filename = figurePNG, 
       width = 297, height = 210, units = "mm")

#### save the figure info ####

dataInfo <- paste(date, "final_report_", fr_name, "_info.TXT", sep = "")

figureInfo <- "3_initiation_zones_windows.R"

source(paste(sourcePath, "4_figure_info.R", sep = ""))

