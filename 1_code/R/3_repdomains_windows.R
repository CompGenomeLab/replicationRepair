#### library ####

library(ggplot2)
library(dplyr)

#### prepare data ####

sourcePath <- paste("~/Documents/My_Projects/Project_Repair_Replication/", 
                    "Scripts/R/", sep = "") 

date <- "[2019.10.31]"

fr_name <- "ERD_LRD_windows"

source(paste(sourcePath, "2_report_sub_dfs.R", sep = ""))

#### set df, variables and rearrange ####

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

#### plot ####

#### filter the data ####

d <- filter(df, phase == "early" | phase == "late", dataset == "ERD")

#### add plot format #### 

source(paste(sourcePath, "4_plot_format.R", sep = ""))

#### create the plot ####

p <- ggplot(d, aes(x = windows, y = log2(xr_ds))) + 
  geom_line(aes(linetype = sample_strand, color = phase)) +
  facet_grid(~product~time_after_exposure~replicate,
             labeller = labeller(product = product_labs, 
                                 time_after_exposure = taex_labs, 
                                 replicate = rep_labs)) +
  xlab(windows_lab) + ylab(fr_xr_ds_lab) +
  scale_x_continuous(limits = c(-half_window-5, half_window+5), 
                     breaks = c(-half_window, 0, half_window), 
                     labels = c("-1 Mbp", "ERDs", "+1 Mbp")) + 
  scale_color_manual(name = "Phase", 
                     label = c("Early Phase", "Late Phase", "Asyncronized"), 
                     values = phase_colors) + 
  labs(color = "Phase", linetype = "Strands")

p <- p + p_format

p # visualize

#### save the plot ####

figurePath <- paste("~/Documents/My_Projects/Project_Repair_Replication/", 
                    "Results/Repdomains/", sep = "")

figureName <- paste(dateout, "Repair_over_Damage_Value_of_", "ERD", 
                    "_for_Every_Sample.pdf", sep = "")

ggsave(path = figurePath, filename = figureName, 
       width = 297, height = 210, units = "mm")

figurePNG <- sub(".pdf", ".png", figureName)

ggsave(path = figurePath, filename = figurePNG, 
       width = 297, height = 210, units = "mm")

#### save the figure info ####

dataInfo <- paste(date, "final_report_", fr_name, "_info.TXT", sep = "")

figureInfo <- "3_repdomains_windows.R"

source(paste(sourcePath, "4_figure_info.R", sep = ""))
