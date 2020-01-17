#### library ####

library(ggplot2)
library(stringr)
library(dplyr)

#### prepare data ####

sourcePath <- paste("~/Documents/My_Projects/Project_Repair_Replication/", 
                    "Scripts/R/", sep = "") 

date <- "[2019.10.31]"

fr_name <- "S1toS2"

source(paste(sourcePath, "2_report_sub_dfs.R", sep = ""))


#### set variables ####

dateout <- Sys.Date()
dateout <- format(dateout, format = "[%Y.%m.%d]")
options(scipen = 999)

df <- fr_xr_ds_ear_la


#### plot ####

#### filter the data ####

#d = filter(df, )

#### add plot format #### 

source(paste(sourcePath, "4_plot_format.R", sep = ""))

#### create the plot ####

p <- ggplot(d) # plot

p <- p +  p_format # adding plot format

p # visualize

#### save the plot ####

figurePath <- paste("~/Documents/My_Projects/Project_Repair_Replication/", 
                    "Results/S1toS2_regions/", sep = "")

figureName <- paste(dateout, "", 
                    fr_name, "", sep = "")

ggsave(path = figurePath, filename = figureName, 
       width = 297, height = 210, units = "mm")

figurePNG <- sub(".pdf", ".png", figureName)

ggsave(path = figurePath, filename = figurePNG, 
       width = 297, height = 210, units = "mm")

#### save the figure info ####

dataInfo <- paste(date, "final_report_", fr_name, "_info.TXT", sep = "")

figureInfo <- "3_S1toS2.R"

source(paste(sourcePath, "4_figure_info.R", sep = ""))


