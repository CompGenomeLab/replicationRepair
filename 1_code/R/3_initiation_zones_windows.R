#### library ####

library(ggplot2)
library(stringr)

#### prepare data ####

sourcePath <- paste("~/Documents/myprojects/replicationRepair/1_code/R/")

date <- "[2019.10.31]"

fr_name <- "InZones_windows"

source(paste(sourcePath, "2_report_sub_dfs.R", sep = ""))

#### set df, variables and rearrange ####

df <- fr

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

d <- filter(df, product == "CPD", method == "Damage_seq")

#### add plot format #### 

source(paste(sourcePath, "4_plot_format.R", sep = ""))

#### create the plot ####

p <- ggplot(d, aes(x = windows, y = RPKM)) + 
  geom_line(aes(linetype = sample_strand, color = sample_strand)) + 
  facet_grid(~product~time_after_exposure~replicate~phase~method, 
             labeller = labeller(product = product_labs, 
                                 method = method_labs, 
                                 time_after_exposure = taex_labs, 
                                 replicate = rep_labs, phase = phase_labs)) + 
  xlab(windows_lab) + ylab(fr_lab) +
  scale_x_continuous(limits = c(-half_window-5, half_window+5), 
                     breaks = c(-half_window, 0, half_window), 
                     labels = c("-10 kb", "Initiation Zones", "+10 kbp")) + 
  scale_color_manual(values = strand_colors) + 
  labs(color = "Strands", linetype = "Strands")

p <- p + p_format

p # visualize

#### save the plot ####

figurePath <- paste("~/Documents/myprojects/replicationRepair/4_output/", 
                    "gitignore/InZones/", sep = "")

figureName <- paste(dateout, "RPKM_Value_of_", fr_name, 
                    "_for_CPD.pdf", sep = "")

ggsave(path = figurePath, filename = figureName, 
       width = 297, height = 210, units = "mm")

figurePNG <- sub(".pdf", ".png", figureName)

ggsave(path = figurePath, filename = figurePNG, 
       width = 297, height = 210, units = "mm")

#### save the figure info ####

dataInfo <- paste(date, "final_report_", fr_name, "_info.TXT", sep = "")

figureInfo <- "3_initiation_zones_windows.R"

source(paste(sourcePath, "4_figure_info.R", sep = ""))


