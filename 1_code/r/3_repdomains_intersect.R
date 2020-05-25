
#### set paths and variables (InZones) ####

# name of the data file and its run date 
fr_name <- "inZones_repdomains_windows_201_100"
date <- "[2020.02.10]"

# path of the codes
sourcePath <- "~/Documents/myprojects/replicationRepair/1_code/R/" 

# save directory of figures
figurePath <- paste("~/Documents/myprojects/replicationRepair/4_output/", 
                    "gitignore/InZones/", sep = "") 

# file name of data information
dataInfo <- paste(date, "final_report_", fr_name, "_info.TXT", sep = "") 

# name of the script
figureInfo <- "3_repdomains_intersect.R" 

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
  df <- domain_name( df, 1 )
  df$dataset <- region
  return(df)
}


#### Repair Rate ####

# rearrange
df <- fr_xr_ds 
df <- rearrange()

phase_labs <- c("Async.", "Early P.", "Late P.")
names(phase_labs) <- c("async", "early", "late")

taex_labs <- c("12 m.", "120 m.")
names(taex_labs) <- c("12", "120")

# filter the data
d <- filter(df, phase == "early" | phase == "late")

# create the plot 
p <- ggplot(d, aes(x = windows, y = log2(xr_ds))) + 
  geom_line(aes(color = sample_strand)) + 
  facet_grid(~product~time_after_exposure~replicate~phase~repdomains, 
             labeller = labeller(product = product_labs, 
                                 time_after_exposure = taex_labs, 
                                 replicate = rep_labs, phase = phase_labs)) + 
  xlab(windows_lab) + ylab(fr_xr_ds_lab) +
  scale_x_continuous(limits = c(-half_window-15, half_window+15), 
                     breaks = c(-half_window, 0, half_window), 
                     labels = c("-10 kb", "Ini. Zones", 
                                "+10 kbp")) + 
  scale_color_manual(values = strand_colors) + 
  labs(color = "Strands")
p <- p + p_format + theme(axis.text.x = element_text(size = 10))
p # visualize

# save figure
figureName <- paste(dateout, "Repair_over_Damage_Value_of_", fr_name, 
                    "_for_Every_Sample.pdf", sep = "")
fig_save( figurePath, figureName )

# save figure info
source(paste(sourcePath, "4_figure_info.R", sep = ""))




#### set paths and variables (SNS-seq) ####

# name of the data file and its run date 
fr_name <- "sns_seq_repdomains_windows_201_100"
date <- "[2020.02.10]"

# path of the codes
sourcePath <- "~/Documents/myprojects/replicationRepair/1_code/R/" 

# save directory of figures
figurePath <- paste("~/Documents/myprojects/replicationRepair/4_output/", 
                    "gitignore/SNS_seq/", sep = "") 

# file name of data information
dataInfo <- paste(date, "final_report_", fr_name, "_info.TXT", sep = "") 

# name of the script
figureInfo <- "3_repdomains_intersect.R" 

# name of the region
region <- "SNS Seq"

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
  df <- domain_name( df, 1 )
  df$dataset <- region
  return(df)
}


#### Repair Rate ####

# rearrange
df <- fr_xr_ds 
df <- rearrange()

phase_labs <- c("Async.", "Early P.", "Late P.")
names(phase_labs) <- c("async", "early", "late")

taex_labs <- c("12 m.", "120 m.")
names(taex_labs) <- c("12", "120")

# filter the data
d <- filter(df, phase == "early" | phase == "late")

# create the plot 
p <- ggplot(d, aes(x = windows, y = log2(xr_ds))) + 
  geom_line(aes(color = sample_strand)) + 
  facet_grid(~product~time_after_exposure~replicate~phase~repdomains, 
             labeller = labeller(product = product_labs, 
                                 time_after_exposure = taex_labs, 
                                 replicate = rep_labs, phase = phase_labs)) + 
  xlab(windows_lab) + ylab(fr_xr_ds_lab) +
  scale_x_continuous(limits = c(-half_window-15, half_window+15), 
                     breaks = c(-half_window, 0, half_window), 
                     labels = c("-10 kb", region,
                                "+10 kbp")) + 
  scale_color_manual(values = strand_colors) + 
  labs(color = "Strands")
p <- p + p_format + theme(axis.text.x = element_text(size = 10))
p # visualize

# save figure
figureName <- paste(dateout, "Repair_over_Damage_Value_of_", fr_name, 
                    "_for_Every_Sample.pdf", sep = "")
fig_save( figurePath, figureName )

# save figure info
source(paste(sourcePath, "4_figure_info.R", sep = ""))




#### set paths and variables (HiRFD) ####

# name of the data file and its run date 
fr_name <- "hiRFD_repdomains_windows_201_100"
date <- "[2020.02.10]"

# path of the codes
sourcePath <- "~/Documents/myprojects/replicationRepair/1_code/R/" 

# save directory of figures
figurePath <- paste("~/Documents/myprojects/replicationRepair/4_output/", 
                    "gitignore/HiRFD/", sep = "") 

# file name of data information
dataInfo <- paste(date, "final_report_", fr_name, "_info.TXT", sep = "") 

# name of the script
figureInfo <- "3_repdomains_intersect.R" 

# name of the region
region <- "High RFDs"

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
  df <- domain_name( df, 1 )
  df$dataset <- region
  return(df)
}


#### Repair Rate ####

# rearrange
df <- fr_xr_ds 
df <- rearrange()

phase_labs <- c("Async.", "Early P.", "Late P.")
names(phase_labs) <- c("async", "early", "late")

taex_labs <- c("12 m.", "120 m.")
names(taex_labs) <- c("12", "120")

# filter the data
d <- filter(df, phase == "early" | phase == "late", repdomains == "ERD")

# create the plot 
p <- ggplot(d, aes(x = windows, y = log2(xr_ds))) + 
  geom_line(aes(color = sample_strand)) + 
  facet_grid(~product~time_after_exposure~replicate~phase~dataset_strand
             ~repdomains, labeller = labeller(product = product_labs, 
                                 time_after_exposure = taex_labs, 
                                 replicate = rep_labs, phase = phase_labs,
                                 dataset_strand = dataset_strand_labs)) + 
  xlab(windows_lab) + ylab(fr_xr_ds_lab) +
  scale_x_continuous(limits = c(-half_window-15, half_window+15), 
                     breaks = c(-half_window, 0, half_window), 
                     labels = c("-10 kb", region, 
                                "+10 kbp")) + 
  scale_color_manual(values = strand_colors) + 
  labs(color = "Strands")
p <- p + p_format + theme(axis.text.x = element_text(size = 10))
p # visualize

# save figure
figureName <- paste(dateout, "Repair_over_Damage_Value_of_", fr_name, 
                    "_for_Every_Sample.pdf", sep = "")
fig_save( figurePath, figureName )

# save figure info
source(paste(sourcePath, "4_figure_info.R", sep = ""))




#### set paths and variables (S1toS2) ####

# name of the data file and its run date 
fr_name <- "S1toS2_repdomains_windows_201_100"
date <- "[2020.02.10]"

# path of the codes
sourcePath <- "~/Documents/myprojects/replicationRepair/1_code/R/" 

# save directory of figures
figurePath <- paste("~/Documents/myprojects/replicationRepair/4_output/", 
                    "gitignore/S1toS2_regions/", sep = "") 

# file name of data information
dataInfo <- paste(date, "final_report_", fr_name, "_info.TXT", sep = "") 

# name of the script
figureInfo <- "3_repdomains_intersect.R" 

# name of the region
region <- "S1toS2"

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
  df <- domain_name( df, 1 )
  df$dataset <- region
  return(df)
}


#### Repair Rate ####

# rearrange
df <- fr_xr_ds 
df <- rearrange()

phase_labs <- c("Async.", "Early P.", "Late P.")
names(phase_labs) <- c("async", "early", "late")

taex_labs <- c("12 m.", "120 m.")
names(taex_labs) <- c("12", "120")

# filter the data
d <- filter(df, phase == "early" | phase == "late")

# create the plot 
p <- ggplot(d, aes(x = windows, y = log2(xr_ds))) + 
  geom_line(aes(color = sample_strand)) + 
  facet_grid(~product~time_after_exposure~replicate~phase~dataset_strand
             ~repdomains, labeller = labeller(product = product_labs, 
                                      time_after_exposure = taex_labs, 
                                      replicate = rep_labs, phase = phase_labs,
                                      dataset_strand = dataset_strand_labs)) +
  xlab(windows_lab) + ylab(fr_xr_ds_lab) +
  scale_x_continuous(limits = c(-half_window-15, half_window+15), 
                     breaks = c(-half_window, 0, half_window), 
                     labels = c("-10 kb", region, 
                                "+10 kbp")) + 
  scale_color_manual(values = strand_colors) + 
  labs(color = "Strands")
p <- p + p_format + theme(axis.text.x = element_text(size = 10))
p # visualize

# save figure
figureName <- paste(dateout, "Repair_over_Damage_Value_of_", fr_name, 
                    "_for_Every_Sample.pdf", sep = "")
fig_save( figurePath, figureName )

# save figure info
source(paste(sourcePath, "4_figure_info.R", sep = ""))



