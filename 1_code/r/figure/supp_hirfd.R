#### set paths and variables ####

repchoice = "B"

# name of the data file and its run date  
fr_name <- "hiRFD_windows_201_10000"
date <- "[2020.02.10]"

# path of the codes
sourcePath <- "~/Documents/myprojects/replicationRepair/1_code/R/" 

# save directory of figures
figurePath <- paste("~/Documents/myprojects/replicationRepair/4_output/", 
                    "gitignore/HiRFD/", sep = "") 

# file name of data information
dataInfo <- paste(date, "final_report_", fr_name, "_info.TXT", sep = "") 

# name of the script
figureInfo <- "figure/supp_hirfd.R" 

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
  df$dataset <- region
  df <- unique(df)
  return(df)
}

method_labs <- c("DS-seq", "XR-seq")
names(method_labs) <- c("Damage_seq", "XR_seq")

#### Repair Rate ####

# rearrange
df <- fr_xr_ds
df <- rearrange()

# filter the data
d <- filter(df, phase != "async", replicate == repchoice)

# create the plot 
p <- p_rr( d )
p <- p + p_format + 
  theme(panel.border = element_rect(fill = NA)) +
  geom_vline(xintercept = 0, color = "gray", linetype = "dashed") +
  facet_grid(~product~time_after_exposure~phase~dataset_strand, 
             labeller = labeller(product = product_labs, 
                                 time_after_exposure = taex_labs, 
                                 replicate = rep_labs, phase = phase_labs,
                                 dataset_strand = dataset_strand_labs))
# p # visualize

# save figure
figureName <- paste(dateout, "Repair_over_Damage_Value_of_", fr_name, 
                    "_rep", repchoice, "_for_Every_Sample.pdf", sep = "")
fig_save( figurePath, figureName )

# save figure info
source(paste(sourcePath, "4_figure_info.R", sep = ""))


#### Early over Late of Repair Rate ####

# rearrange
df <- fr_xr_ds_ear_la
df <- rearrange()

# filter the data
d <- filter(df, replicate == repchoice)

# create the plot 
p <- p_rr_el( d )
p <- p + p_format +
  theme(panel.border = element_rect(fill = NA)) +
  geom_vline(xintercept = 0, color = "gray", linetype = "dashed") +
  facet_grid(~product~time_after_exposure~dataset_strand,
             labeller = labeller(product = product_labs, 
                                 time_after_exposure = taex_labs, 
                                 replicate = rep_labs,
                                 dataset_strand = dataset_strand_labs))
# p # visualize

# save figure
figureName <- paste(dateout, "Early_over_Late_Phase_Value_of_", 
                    "Repair_over_Damage_Value_of_", fr_name, 
                    "_rep", repchoice, "_for_Every_Sample.pdf", sep = "")
fig_save( figurePath, figureName )

# save figure info
source(paste(sourcePath, "4_figure_info.R", sep = ""))



#### Plus over Minus of Repair Rate ####

# rearrange
df <- fr_xr_ds_plus_min
df <- rearrange()

# filter the data
d <- filter(df, phase != "async", replicate == repchoice)

# create the plot 
p <- p_rr_pm( d )
p <- p + p_format +
  theme(panel.border = element_rect(fill = NA)) +
  geom_vline(xintercept = 0, color = "gray", linetype = "dashed") +
  facet_grid(~product~time_after_exposure~dataset_strand, 
             labeller = labeller(product = product_labs, 
                                 time_after_exposure = taex_labs, 
                                 replicate = rep_labs, 
                                 dataset_strand = dataset_strand_labs))
# p # visualize

# save figure
figureName <- paste(dateout, "Plus_over_Minus_Strand_Value_of_", 
                    "Repair_over_Damage_Value_of_", fr_name, 
                    "_rep", repchoice, "_for_Every_Sample.pdf", sep = "")
fig_save( figurePath, figureName )

# save figure info
source(paste(sourcePath, "4_figure_info.R", sep = ""))


#### Early over Late ####

# rearrange
df <- fr_ear_la
df <- rearrange()

# filter the data
d <- filter(df, method != "DNA_seq", replicate == repchoice )

# create the plot 
p <- p_el( d )
p <- p + p_format +
  theme(panel.border = element_rect(fill = NA)) +
  geom_vline(xintercept = 0, color = "gray", linetype = "dashed") +
  facet_grid(~product~time_after_exposure~method~dataset_strand,
             labeller = labeller(product = product_labs, 
                                 time_after_exposure = taex_labs, 
                                 replicate = rep_labs,
                                 method = method_labs, 
                                 dataset_strand = dataset_strand_labs)) 
# p # visualize

# save figure
figureName <- paste(dateout, "Early_over_Late_Phase_Value_of_", 
                    fr_name, "_rep", repchoice, "_for_Every_Sample.pdf", sep = "")
fig_save( figurePath, figureName )

# save figure info
source(paste(sourcePath, "4_figure_info.R", sep = ""))


#### Plus over Minus ####

# rearrange
df <- fr_plus_min
df <- rearrange()

# filter the data
d <- filter(df, method != "DNA_seq", phase != "async", replicate == repchoice)

# create the plot 
p <- p_pm( d )
p <- p + p_format +
  theme(panel.border = element_rect(fill = NA)) +
  geom_vline(xintercept = 0, color = "gray", linetype = "dashed") +
  facet_grid(~product~time_after_exposure~method~dataset_strand, 
             labeller = labeller(product = product_labs, 
                                 method = method_labs, 
                                 time_after_exposure = taex_labs, 
                                 replicate = rep_labs,
                                 dataset_strand = dataset_strand_labs))
# p # visualize

# save figure
figureName <- paste(dateout, "Plus_over_Minus_Value_of_", fr_name, 
                    "_rep", repchoice, "_for_Every_Sample.pdf", sep = "")
fig_save( figurePath, figureName )

# save figure info
source(paste(sourcePath, "4_figure_info.R", sep = ""))


