
library(patchwork)

#### set paths and variables ####

# name of the data file and its run date  
fr_name <- "inZones_windows_201_100"
date <- "[2020.02.10]"

# path of the codes
sourcePath <- "~/Documents/myprojects/replicationRepair/1_code/R/" 

# save directory of figures
figurePath <- paste("~/Documents/myprojects/replicationRepair/4_output/", 
                    "gitignore/InZones/", sep = "") 

# file name of data information
dataInfo <- paste(date, "final_report_", fr_name, "_info.TXT", sep = "") 

# name of the script
figureInfo <- "3_initiation_zones_windows.R" 

# name of the region
region <- "Initiation Zones"

# pre-analysis of data file
source(paste(sourcePath, "2_report_sub_dfs.R", sep = ""))

# functions
source(paste(sourcePath, "4_functions.R", sep = ""))

# add plot format
source(paste(sourcePath, "4_plot_format.R", sep = ""))

# date of today
dateout <- output_date()

# all the functions that will be used for rearrangement
rearrange <- function( ){
  window_number <- as.numeric(get_sub_str( fr_name, "windows_", "_1" ))
  window_length <- as.numeric(get_sub_str( fr_name, paste(window_number, 
                                                          "_", sep = ""), "$" ))
  half_window <<- middle( window_number ) 
  rlength <<- region_length( window_number, window_length )
  df <- window_numbering( df, 4, half_window )
  df$dataset <- region
  df <- unique(df)
  return(df)
}


#### p.A.1.1 and p.A.2.1 ####

# rearrange
df <- fr
df <- rearrange()

# filter the data
d <- filter(df, method != "DNA_seq", phase != "async", replicate == "B", 
            method == "Damage_seq")

# create the plot 
p.A.1.1 <- p_RPKM( d )
p.A.1.1 <- p.A.1.1 + p_format + 
  geom_vline(xintercept = 0, color = "gray", linetype = "dashed") +
  ylim(0.10, 0.25) +
  facet_grid(~product~time_after_exposure~phase~method, 
             labeller = labeller(product = product_labs, 
                                 method = method_labs, 
                                 time_after_exposure = taex_labs,
                                 phase = phase_labs)) + 
  scale_x_continuous(limits = c(-half_window, half_window), 
                     breaks = c(-half_window, 0, half_window), 
                     labels = c("-10", "0", "+10")) + 
  xlab("") +
  theme(strip.text.y = element_blank(),
        panel.border = element_rect(fill = NA)) 

p.A.1.1 # visualize

# filter the data
d <- filter(df, method != "DNA_seq", phase != "async", replicate == "B", 
            method == "XR_seq")

# create the plot 
p.A.2.1 <- p_RPKM( d )
p.A.2.1 <- p.A.2.1 + p_format + 
  geom_vline(xintercept = 0, color = "gray", linetype = "dashed") +
  ylim(0.10, 0.35) +
  facet_grid(~product~time_after_exposure~phase~method, 
             labeller = labeller(product = product_labs, 
                                 method = method_labs, 
                                 time_after_exposure = taex_labs,
                                 phase = phase_labs)) + 
  scale_x_continuous(limits = c(-half_window, half_window), 
                     breaks = c(-half_window, 0, half_window), 
                     labels = c("-10", "0", "+10")) + 
  xlab("") +
  theme(strip.text.y = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_rect(fill = NA)) 

p.A.2.1 # visualize



#### p.A.2.2 xr sim ####

#### library 

library(ggplot2)
library(stringr)

#### prepare data

sourcePath <- paste("~/Documents/myprojects/replicationRepair/1_code/R/")

date <- "[2020.02.25]"

fr_name <- "inZones_windows_201_100_sim"

source(paste(sourcePath, "2_report_sub_dfs.R", sep = ""))

#### set df, variables and rearrange 

df <- fr

for (rearrange in 1) {  
  
  #### set variables 
  
  dateout <- Sys.Date()
  dateout <- format(dateout, format = "[%Y.%m.%d]")
  
  window_number <- 201
  if (window_number %% 2 == 0) {
    half_window <- window_number / 2
  } else {
    half_window <- (window_number - 1) / 2 + 1
  }
  
  #### rearrange 
  
  # window numbers separated from dataset names
  window <- data.frame(str_split_fixed(df$dataset, "_", -1))
  
  df$windows <- as.numeric(levels(
    window[ , ncol(window)]))[window[ , ncol(window)]] - half_window
  
  df$dataset <- "Initiation Zones"
  
  rm(window, rearrange)
  
  
}

#### plot 

#### filter the data 

d <- filter(df, phase != "async", method != "DNA_seq", replicate == "B")

#### add plot format  

source(paste(sourcePath, "4_plot_format.R", sep = ""))

method_labs <- c("Simulated \nXR-seq")
names(method_labs) <- c("XR_seq")

#### create the plot 

p.A.2.2 <- ggplot(d, aes(x = windows, y = RPKM)) + 
  geom_line(aes(color = sample_strand)) + 
  geom_vline(xintercept = 0, color = "gray", linetype = "dashed") +
  facet_grid(~product~time_after_exposure~phase~method, 
             labeller = labeller(product = product_labs, 
                                 method = method_labs, 
                                 time_after_exposure = taex_labs, 
                                 phase = phase_labs)) + 
  xlab("Relative Position Centered at Initiation Zones (kb)") +
  ylab(fr_lab) +
  scale_x_continuous(limits = c(-half_window, half_window), 
                     breaks = c(-half_window, 0, half_window), 
                     labels = c("-10", "0", "+10")) + 
  scale_color_manual(values = strand_colors) + 
  labs(color = "Strands")

p.A.2.2 <- p.A.2.2 + p_format + 
  ylim(0.10, 0.35) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.border = element_rect(fill = NA)
  )
p.A.2.2 # visualize

#### p.A.1.2 ds sim ####

#### library

library(ggplot2)
library(stringr)

#### prepare data

sourcePath <- paste("~/Documents/myprojects/replicationRepair/1_code/R/")

date <- "[2020.01.19]"

fr_name <- "InZones_windows_sim"

source(paste(sourcePath, "2_report_sub_dfs.R", sep = ""))

#### set df, variables and rearrange 

df <- fr

for (rearrange in 1) {  
  
  #### set variables
  
  dateout <- Sys.Date()
  dateout <- format(dateout, format = "[%Y.%m.%d]")
  
  window_number <- 201
  if (window_number %% 2 == 0) {
    half_window <- window_number / 2
  } else {
    half_window <- (window_number - 1) / 2 + 1
  }
  
  #### rearrange
  
  # window numbers separated from dataset names
  window <- data.frame(str_split_fixed(df$dataset, "_", -1))
  
  df$windows <- as.numeric(levels(
    window[ , ncol(window)]))[window[ , ncol(window)]] - half_window
  
  df$dataset <- "Initiation Zones"
  
  rm(window, rearrange)
  
  
}

#### plot 

#### filter the data 

d <- filter(df, phase != "async", replicate == "B")

#### add plot format 

source(paste(sourcePath, "4_plot_format.R", sep = ""))

method_labs <- c("Simulated \nDamage-seq")
names(method_labs) <- c("Damage_seq")

#### create the plot 

p.A.1.2 <- ggplot(d, aes(x = windows, y = RPKM)) + 
  geom_line(aes(color = sample_strand)) + 
  geom_vline(xintercept = 0, color = "gray", linetype = "dashed") +
  facet_grid(~product~time_after_exposure~phase~method,
             labeller = labeller(product = product_labs,
                                 method = method_labs,
                                 time_after_exposure = taex_labs,
                                 phase = phase_labs)) +
  xlab("Relative Position Centered at Initiation Zones (kb)") + 
  ylab(fr_lab) +
  scale_x_continuous(limits = c(-half_window, half_window), 
                     breaks = c(-half_window, 0, half_window), 
                     labels = c("-10", "0", "+10")) + 
  scale_color_manual(values = strand_colors) + 
  labs(color = "Strands")

p.A.1.2 <- p.A.1.2 + p_format + 
  ylim(0.10, 0.25) +
  theme(axis.title.y=element_blank(),
        panel.border = element_rect(fill = NA),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
         strip.text.y = element_blank())

p.A.1.2 # visualize

#### final plot ####

layout <- "
AABBCCDD
"

p.A.1.1 + p.A.1.2 + p.A.2.1 + p.A.2.2 +
  plot_layout(design = layout) +
  #plot_annotation(tag_levels = c('A'), tag_suffix = ':') +
  plot_layout(guides = "collect", tag_level = 'new') & 
  theme(legend.position = 'bottom') 

ggsave("~/Desktop/fig4_repB.png", width = 297, height = 210, units = "mm")







