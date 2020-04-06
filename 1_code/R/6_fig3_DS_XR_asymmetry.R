
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
  return(df)
}


#### p.A.1.1 and p.A.2.1 ####

# rearrange
df <- fr
df <- rearrange()

# filter the data
d <- filter(df, method != "DNA_seq", phase != "async", replicate == "A", 
            method == "Damage_seq")

# create the plot 
p.A.1.1 <- p_RPKM( d )
p.A.1.1 <- p.A.1.1 + p_format + 
  facet_grid(~product~time_after_exposure~phase~method, 
             labeller = labeller(product = product_labs, 
                                 method = method_labs, 
                                 time_after_exposure = taex_labs,
                                 phase = phase_labs)) + 
  ylim(0.10, 0.25) + 
  theme(strip.text.y = element_blank())

p.A.1.1 # visualize

 # filter the data
 d <- filter(df, method != "DNA_seq", phase != "async", replicate == "A", 
             method == "XR_seq")
 
 # create the plot 
 p.A.2.1 <- p_RPKM( d )
 p.A.2.1 <- p.A.2.1 + p_format + 
   facet_grid(~product~time_after_exposure~phase~method, 
              labeller = labeller(product = product_labs, 
                                  method = method_labs, 
                                  time_after_exposure = taex_labs,
                                  phase = phase_labs)) + 
   ylim(0.10, 0.35) + 
   theme(strip.text.y = element_blank())

  p.A.2.1 # visualize


#### Plus over Minus ####

# rearrange
df <- fr_plus_min
df <- rearrange()

# filter the data
d <- filter(df, phase != "async", method != "DNA_seq", replicate == "A")

# create the plot 
p <- p_mp( d )
p.C <- p + p_format + 
  facet_grid(~product~time_after_exposure~method, 
              labeller = labeller(product = product_labs, 
                                  method = method_labs, 
                                  time_after_exposure = taex_labs))


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

d <- filter(df, phase != "async", method != "DNA_seq", replicate == "A")

#### add plot format  

source(paste(sourcePath, "4_plot_format.R", sep = ""))

method_labs <- c("Simulated XR-seq")
names(method_labs) <- c("XR_seq")

#### create the plot 

p.A.2.2 <- ggplot(d, aes(x = windows, y = RPKM)) + 
  geom_line(aes(color = sample_strand)) + 
  facet_grid(~product~time_after_exposure~phase~method, 
             labeller = labeller(product = product_labs, 
                                 method = method_labs, 
                                 time_after_exposure = taex_labs, 
                                 phase = phase_labs)) + 
  xlab(windows_lab) + ylab(fr_lab) +
  scale_x_continuous(limits = c(-half_window, half_window), 
                     breaks = c(-half_window, 0, half_window), 
                     labels = c("-10 kb", "0", "+10 kbp")) + 
  scale_color_manual(values = strand_colors) + 
  labs(color = "Strands")

p.A.2.2 <- p.A.2.2 + p_format + ylim(0.10, 0.35) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

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

d <- filter(df, phase != "async", replicate == "A")

#### add plot format 

source(paste(sourcePath, "4_plot_format.R", sep = ""))

method_labs <- c("Simulated Damage-seq")
names(method_labs) <- c("Damage_seq")

#### create the plot 

p.A.1.2 <- ggplot(d, aes(x = windows, y = RPKM)) + 
  geom_line(aes(color = sample_strand)) + 
  facet_grid(~product~time_after_exposure~phase~method, 
             labeller = labeller(product = product_labs, 
                                 method = method_labs, 
                                 time_after_exposure = taex_labs, 
                                 phase = phase_labs)) + 
  xlab(windows_lab) + ylab(fr_lab) +
  scale_x_continuous(limits = c(-half_window, half_window), 
                     breaks = c(-half_window, 0, half_window), 
                     labels = c("-10 kb", "0", "+10 kbp")) + 
  scale_color_manual(values = strand_colors) + 
  labs(color = "Strands")

p.A.1.2 <- p.A.1.2 + p_format + ylim(0.10, 0.25) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        strip.text.y = element_blank())

p.A.1.2 # visualize

#### #####

Ccontent <- read.table(
  paste(mut_dir, 
        "[2020.03.04]inZones_nuc_content_norepdomain_201_100.txt", 
        sep = ""), header = TRUE)

Ccontent <- Ccontent[,c(1,6,7)]

names(Ccontent) <- c("name","+", "-")

Ccontent_plus <- Ccontent[,c(1,2)]
Ccontent_minus <- Ccontent[,c(1,3)]

Ccontent_plus <- melt(Ccontent_plus, measure.vars = "+")
Ccontent_minus <- melt(Ccontent_minus, measure.vars = "-")

Ccontent <- rbind(Ccontent_plus, Ccontent_minus)

names(Ccontent) <- c("name","strands", "C_content")

for (rearrange in 1) {  
  
  window_number <- 201
  if (window_number %% 2 == 0) {
    half_window <- window_number / 2
  } else {
    half_window <- (window_number - 1) / 2 + 1
  }
  
  #### rearrange
  
  # window numbers separated from dataset names
  window <- data.frame(str_split_fixed(Ccontent$name, "_", -1))
  
  Ccontent$windows <- as.numeric(levels(
    window[ , ncol(window)]))[window[ , ncol(window)]] - half_window
  
  Ccontent$name <- "Initiation Zones"
  
  rm(window, rearrange)
  
  
}

p.B.1.3 <- ggplot(Ccontent, aes(x = windows, y = as.numeric(C_content))) + 
  geom_line(aes( color = strands )) + p_format +
  xlab(windows_lab) + ylab("C Content Frequency (%) \nat Initiation Zones") +
  scale_x_continuous(limits = c(-100, 100), 
                     breaks = c(-100, 0, 100), 
                     labels = c(paste("-10 kb", sep = ""), 
                                "0", 
                                paste("+10 kb", sep = ""))) + 
  scale_color_manual(values = strand_colors, guide = FALSE) 
p.B.1.3

Ccontent_plus <- filter(Ccontent, windows > 0)
Ccontent_plus_agg <- aggregate(x = Ccontent_plus$C_content, by = 
                            list(Ccontent_plus$name, Ccontent_plus$strands), 
                            FUN = "mean")
Ccontent_plus_agg$direction <- "Right Replicating"
Ccontent_minus <- filter(Ccontent, windows < 0)
Ccontent_minus_agg <- aggregate(x = Ccontent_minus$C_content, by = 
                             list(Ccontent_minus$name, Ccontent_minus$strands), 
                             FUN = "mean")
Ccontent_minus_agg$direction <- "Left Replicating"
Ccontent_agg <- rbind(Ccontent_plus_agg, Ccontent_minus_agg)


p.B.1.1 <- ggplot() + 
  geom_bar(data = Ccontent_agg, aes(x = Group.2, y = x, 
                               fill = Group.2), 
           stat = "identity", size = 1.5, 
           position=position_dodge2(padding = 0.05)) +
  facet_wrap(~direction) +
  xlab("Replication Domains") + 
  ylab("Cytosine Count per Region") +
  scale_fill_manual(values = c("+" = "#0571b0", 
                                "-" = "#ca0020"), 
                     guide = FALSE) +
  labs(color = "Strands", fill = "") 

mut_casted <- dcast(Ccontent, name + windows ~ strands, 
                    value.var = "C_content")
mut_casted$plus_min <- mut_casted$"+" / mut_casted$"-" 
mut_plus <- filter(mut_casted, windows > 0)
mut_plus_agg <- aggregate(x = mut_plus$plus_min, by = 
                            list(mut_plus$name), 
                          FUN = "mean")
mut_plus_agg$direction <- "Right Replicating"
mut_minus <- filter(mut_casted, windows < 0)
mut_minus_agg <- aggregate(x = mut_minus$plus_min, by = 
                             list(mut_minus$name), 
                           FUN = "mean")
mut_minus_agg$direction <- "Left Replicating"
mut_agg <- rbind(mut_plus_agg, mut_minus_agg)

#### create the plot

p.B.1.2 <- ggplot() + 
  geom_bar(data = mut_agg, aes(x = direction, y = log2(x), fill = direction), 
           stat = "identity", position=position_dodge()) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  xlab("Replication Direction") + ylab("Plus/Minus \n(Log2)") +
  scale_fill_manual(values = c("gray", "gray"), guide = FALSE) 

# add plot format
p.B.1.2 <- p.B.1.2 + p_format + 
  theme(strip.background = element_blank(),
        strip.text.x = element_blank())  

p.B.1.2 # visualize

#### final plot ####

layout <- "
AABBCCDD
AABBCCDD
FFFFGGGG
#EE#GGGG
"

(p.A.1.1 + p.A.1.2 + p.A.2.1) + p.A.2.2 + p.B.1.2 + p.B.1.3 + p.C +
  plot_layout(design = layout) 
  #plot_annotation(tag_levels = c('A'), tag_suffix = ':')

ggsave("~/Desktop/fig3.png", width = 497, height = 410, units = "mm")







