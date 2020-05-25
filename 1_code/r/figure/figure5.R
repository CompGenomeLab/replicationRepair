library(patchwork)

#### set ####

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


#### Repair Rate ####

# rearrange
df <- fr_xr_ds
df <- rearrange()

#### p.A.1 ####

# filter the data
d <- filter(df, replicate == "B", 
            product == "64_PP", time_after_exposure == "12", phase == "early")

# create the plot 
p <- p_rr( d )
p.A.1 <- p + p_format + ylim(-0.1, 0.7) +
  scale_x_continuous(limits = c(-half_window, half_window), 
                     breaks = c(-half_window, 0, half_window), 
                     labels = c("-10", "0", "+10")) + 
  geom_vline(xintercept = 0, color = "gray", linetype = "dashed") +
  facet_grid(~replicate) + theme(strip.background = element_blank(),
                                 strip.text = element_blank(),
                                 axis.title.x=element_blank())

p.A.1 # visualize

#### p.B.1 ####

# filter the data
d <- filter(df, phase != "async", replicate == "B", 
            product == "64_PP", time_after_exposure == "12", phase == "late")

# create the plot 
p <- p_rr( d )
p.B.1 <- p + p_format + ylim(-0.1, 0.7) +
  scale_x_continuous(limits = c(-half_window, half_window), 
                     breaks = c(-half_window, 0, half_window), 
                     labels = c("-10", "0", "+10")) + 
  geom_vline(xintercept = 0, color = "gray", linetype = "dashed") +
  facet_grid(~replicate) + theme(strip.background = element_blank(),
                                 strip.text = element_blank(),
                                 axis.title.x=element_blank())

p.B.1 # visualize

#### p.C.1 ####

# filter the data
d <- filter(df, phase != "async", replicate == "B", 
            product == "CPD", time_after_exposure == "12", phase == "early")

# create the plot 
p <- p_rr( d )
p.C.1 <- p + p_format + ylim(-0.1, 0.7) +
  scale_x_continuous(limits = c(-half_window, half_window), 
                     breaks = c(-half_window, 0, half_window), 
                     labels = c("-10", "0", "+10")) + 
  geom_vline(xintercept = 0, color = "gray", linetype = "dashed") +
  facet_grid(~replicate) + theme(strip.background = element_blank(),
                                 strip.text = element_blank(),
                                 axis.title.x=element_blank())

p.C.1 # visualize

#### p.D.1 ####

# filter the data
d <- filter(df, phase != "async", replicate == "B", 
            product == "CPD", time_after_exposure == "12", phase == "late")

# create the plot 
p <- p_rr( d )
p.D.1 <- p + p_format + ylim(-0.1, 0.7) +
  scale_x_continuous(limits = c(-half_window, half_window), 
                     breaks = c(-half_window, 0, half_window), 
                     labels = c("-10", "0", "+10")) + 
  geom_vline(xintercept = 0, color = "gray", linetype = "dashed") +
  facet_grid(~replicate) + theme(strip.background = element_blank(),
                                 strip.text = element_blank(),
                                 axis.title.x=element_blank())

p.D.1 # visualize

#### p.E.1 ####

# filter the data
d <- filter(df, phase != "async", replicate == "B", 
            product == "CPD", time_after_exposure == "120", phase == "early")

# create the plot 
p <- p_rr( d )
p.E.1 <- p + p_format + ylim(-0.1, 0.7) +
  scale_x_continuous(limits = c(-half_window, half_window), 
                     breaks = c(-half_window, 0, half_window), 
                     labels = c("-10", "0", "+10")) + 
  geom_vline(xintercept = 0, color = "gray", linetype = "dashed") +
  facet_grid(~replicate) + theme(strip.background = element_blank(),
                                 strip.text = element_blank(),
                                 axis.title.x=element_blank())

p.E.1 # visualize

#### p.F.1 ####

# filter the data
d <- filter(df, phase != "async", replicate == "B", 
            product == "CPD", time_after_exposure == "120", phase == "late")

# create the plot 
p <- p_rr( d )
p.F.1 <- p + p_format + ylim(-0.1, 0.7) +
  scale_x_continuous(limits = c(-half_window, half_window), 
                     breaks = c(-half_window, 0, half_window), 
                     labels = c("-10", "0", "+10")) + 
  xlab("Relative Position Centered \nat Initiation Zones (kb)") + 
  geom_vline(xintercept = 0, color = "gray", linetype = "dashed") +
  facet_grid(~replicate) + theme(strip.background = element_blank(),
                                 strip.text = element_blank())

p.F.1 # visualize

#### set ####

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
  df <- unique(df)
  return(df)
}


#### Repdomain Repair Rate ####

df <- fr_xr_ds 
df <- rearrange()

#### p.A.2 ####

phase_labs <- c("Async.", "Early P.", "Late P.")
names(phase_labs) <- c("async", "early", "late")

taex_labs <- c("12 m.", "120 m.")
names(taex_labs) <- c("12", "120")

# filter the data
d <- filter(df, phase != "async", replicate == "B", 
            product == "64_PP", time_after_exposure == "12", phase == "early")

# create the plot 
p <- ggplot(d, aes(x = windows, y = log2(xr_ds))) + 
  geom_line(aes(color = sample_strand)) +  
  facet_grid(~repdomains) +
  xlab(windows_lab) + ylab(fr_xr_ds_lab) +
  scale_x_continuous(limits = c(-half_window, half_window), 
                     breaks = c(-half_window, 0, half_window), 
                     labels = c("-10", "0", "+10")) + 
  scale_color_manual(values = strand_colors) + 
  labs(color = "Strands")
p.A.2 <- p + p_format + ylim(-1.5, 1.5) + 
  geom_vline(xintercept = 0, color = "gray", linetype = "dashed") +
  theme(axis.title.x=element_blank(),
                                                axis.title.y=element_blank(),
                                                #axis.text.y=element_blank(),
                                                axis.text.x=element_text(hjust=c(0.3, 0.5, 0.7)))

p.A.2 # visualize


#### p.A.3 ####

mut_plus <- filter(d, windows > 0)
mut_plus_agg <- aggregate(x = mut_plus$xr_ds, by = 
                            list(mut_plus$dataset, mut_plus$repdomains, 
                                 mut_plus$sample_strand), FUN = "mean")
mut_plus_agg$direction <- "Right Replicating"
mut_minus <- filter(d, windows < 0)
mut_minus_agg <- aggregate(x = mut_minus$xr_ds, by = 
                             list(mut_minus$dataset, mut_minus$repdomains, 
                                  mut_minus$sample_strand), FUN = "mean")
mut_minus_agg$direction <- "Left Replicating"
mut_agg <- rbind(mut_plus_agg, mut_minus_agg)

#### create the plot 

p.A.3 <- ggplot() + 
  geom_bar(data = mut_agg, aes(x = Group.2, y = x, 
                               fill = Group.3), 
           stat = "identity", size = 1.5, 
           position=position_dodge2(padding = 0.05)) +
  facet_wrap(~direction) +
  xlab("Replication Domains") + 
  ylab("Repair \nRate") +
  scale_fill_manual(values = c("+" = "#0571b0", 
                               "-" = "#ca0020"), 
                    guide = FALSE) +
  labs(color = "Strands", fill = "") 

p.A.3 <- p.A.3 + p_format + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) + ylim(0, 2.5) # add plot format

p.A.3 # visualize

#### p.A.4 ####

mut_casted <- dcast(d, dataset + repdomains + windows ~ sample_strand, 
                    value.var = "xr_ds")
mut_casted$plus_min <- mut_casted$"+" / mut_casted$"-" 
mut_plus <- filter(mut_casted, windows > 0)
mut_plus_agg <- aggregate(x = mut_plus$plus_min, by = 
                            list(mut_plus$dataset, mut_plus$repdomains), 
                          FUN = "mean")
mut_plus_agg$direction <- "Right Replicating"
mut_minus <- filter(mut_casted, windows < 0)
mut_minus_agg <- aggregate(x = mut_minus$plus_min, by = 
                             list(mut_minus$dataset, mut_minus$repdomains), 
                           FUN = "mean")
mut_minus_agg$direction <- "Left Replicating"
mut_agg <- rbind(mut_plus_agg, mut_minus_agg)

#### create the plot

p.A.4 <- ggplot() + 
  geom_bar(data = mut_agg, aes(x = Group.2, y = log2(x)), 
           stat = "identity", position=position_dodge()) +
  facet_wrap(~direction) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  xlab("") + ylab("Plus/Minus \n(Log2)") +
  scale_fill_manual(values = repdomain_colors, guide = FALSE) 

# add plot format
p.A.4 <- p.A.4 + p_format + 
  theme(strip.background = element_blank(),
        strip.text.x = element_blank()) + ylim(-0.4, 0.4) 

p.A.4 # visualize


#### p.B.2 ####

# filter the data
d <- filter(df, phase != "async", replicate == "B", 
            product == "64_PP", time_after_exposure == "12", phase == "late")

# create the plot 
p <- ggplot(d, aes(x = windows, y = log2(xr_ds))) + 
  geom_line(aes(color = sample_strand)) + 
  facet_grid(~repdomains) +
  xlab(windows_lab) + ylab(fr_xr_ds_lab) +
  scale_x_continuous(limits = c(-half_window, half_window), 
                     breaks = c(-half_window, 0, half_window), 
                     labels = c("-10", "0", "+10")) + 
  scale_color_manual(values = strand_colors) + 
  labs(color = "Strands")
p.B.2 <- p + p_format + ylim(-1.5, 1.5) + 
  geom_vline(xintercept = 0, color = "gray", linetype = "dashed") +
  theme(axis.title.x=element_blank(),
                                                axis.title.y=element_blank(),
                                                #axis.text.y=element_blank(),
                                                axis.text.x=element_text(hjust=c(0.3, 0.5, 0.7)))

p.B.2 # visualize



#### p.B.3 ####

mut_plus <- filter(d, windows > 0)
mut_plus_agg <- aggregate(x = mut_plus$xr_ds, by = 
                            list(mut_plus$dataset, mut_plus$repdomains, 
                                 mut_plus$sample_strand), FUN = "mean")
mut_plus_agg$direction <- "Right Replicating"
mut_minus <- filter(d, windows < 0)
mut_minus_agg <- aggregate(x = mut_minus$xr_ds, by = 
                             list(mut_minus$dataset, mut_minus$repdomains, 
                                  mut_minus$sample_strand), FUN = "mean")
mut_minus_agg$direction <- "Left Replicating"
mut_agg <- rbind(mut_plus_agg, mut_minus_agg)

#### create the plot 

p.B.3 <- ggplot() + 
  geom_bar(data = mut_agg, aes(x = Group.2, y = x, 
                               fill = Group.3), 
           stat = "identity", size = 1.5, 
           position=position_dodge2(padding = 0.05)) +
  facet_wrap(~direction) +
  xlab("Replication Domains") + 
  ylab("Repair \nRate") +
  scale_fill_manual(values = c("+" = "#0571b0", 
                               "-" = "#ca0020"), 
                    guide = FALSE) +
  labs(color = "Strands", fill = "") 

p.B.3 <- p.B.3 + p_format + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) + ylim(0, 2.5) # add plot format

#### p.B.4 ####

mut_casted <- dcast(d, dataset + repdomains + windows ~ sample_strand, 
                    value.var = "xr_ds")
mut_casted$plus_min <- mut_casted$"+" / mut_casted$"-" 
mut_plus <- filter(mut_casted, windows > 0)
mut_plus_agg <- aggregate(x = mut_plus$plus_min, by = 
                            list(mut_plus$dataset, mut_plus$repdomains), 
                          FUN = "mean")
mut_plus_agg$direction <- "Right Replicating"
mut_minus <- filter(mut_casted, windows < 0)
mut_minus_agg <- aggregate(x = mut_minus$plus_min, by = 
                             list(mut_minus$dataset, mut_minus$repdomains), 
                           FUN = "mean")
mut_minus_agg$direction <- "Left Replicating"
mut_agg <- rbind(mut_plus_agg, mut_minus_agg)

#### create the plot

p.B.4 <- ggplot() + 
  geom_bar(data = mut_agg, aes(x = Group.2, y = log2(x)), 
           stat = "identity", position=position_dodge()) +
  facet_wrap(~direction) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  xlab("") + ylab("Plus/Minus \n(Log2)") +
  scale_fill_manual(values = repdomain_colors, guide = FALSE) 

# add plot format
p.B.4 <- p.B.4 + p_format + 
  theme(strip.background = element_blank(),
        strip.text.x = element_blank()) + ylim(-0.4, 0.4) 

p.B.4 # visualize


#### p.C.2 ####

# filter the data
d <- filter(df, phase != "async", replicate == "B", 
            product == "CPD", time_after_exposure == "12", phase == "early")

# create the plot 
p <- ggplot(d, aes(x = windows, y = log2(xr_ds))) + 
  geom_line(aes(color = sample_strand)) + 
  facet_grid(~repdomains) +
  xlab(windows_lab) + ylab(fr_xr_ds_lab) +
  scale_x_continuous(limits = c(-half_window, half_window), 
                     breaks = c(-half_window, 0, half_window), 
                     labels = c("-10", "0", "+10")) + 
  scale_color_manual(values = strand_colors) + 
  labs(color = "Strands")
p.C.2 <- p + p_format + ylim(-1.5, 1.5) + 
  geom_vline(xintercept = 0, color = "gray", linetype = "dashed") +
  theme(axis.title.x=element_blank(),
                                                axis.title.y=element_blank(),
                                                #axis.text.y=element_blank(),
                                                axis.text.x=element_text(hjust=c(0.3, 0.5, 0.7)))

p.C.2 # visualize



#### p.C.3 ####

mut_plus <- filter(d, windows > 0)
mut_plus_agg <- aggregate(x = mut_plus$xr_ds, by = 
                            list(mut_plus$dataset, mut_plus$repdomains, 
                                 mut_plus$sample_strand), FUN = "mean")
mut_plus_agg$direction <- "Right Replicating"
mut_minus <- filter(d, windows < 0)
mut_minus_agg <- aggregate(x = mut_minus$xr_ds, by = 
                             list(mut_minus$dataset, mut_minus$repdomains, 
                                  mut_minus$sample_strand), FUN = "mean")
mut_minus_agg$direction <- "Left Replicating"
mut_agg <- rbind(mut_plus_agg, mut_minus_agg)

#### create the plot 

p.C.3 <- ggplot() + 
  geom_bar(data = mut_agg, aes(x = Group.2, y = x, 
                               fill = Group.3), 
           stat = "identity", size = 1.5, 
           position=position_dodge2(padding = 0.05)) +
  facet_wrap(~direction) +
  xlab("Replication Domains") + 
  ylab("Repair \nRate") +
  scale_fill_manual(values = c("+" = "#0571b0", 
                               "-" = "#ca0020"), 
                    guide = FALSE) +
  labs(color = "Strands", fill = "") 

p.C.3 <- p.C.3 + p_format + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) + ylim(0, 2.5) # add plot format

#### p.C.4 ####

mut_casted <- dcast(d, dataset + repdomains + windows ~ sample_strand, 
                    value.var = "xr_ds")
mut_casted$plus_min <- mut_casted$"+" / mut_casted$"-" 
mut_plus <- filter(mut_casted, windows > 0)
mut_plus_agg <- aggregate(x = mut_plus$plus_min, by = 
                            list(mut_plus$dataset, mut_plus$repdomains), 
                          FUN = "mean")
mut_plus_agg$direction <- "Right Replicating"
mut_minus <- filter(mut_casted, windows < 0)
mut_minus_agg <- aggregate(x = mut_minus$plus_min, by = 
                             list(mut_minus$dataset, mut_minus$repdomains), 
                           FUN = "mean")
mut_minus_agg$direction <- "Left Replicating"
mut_agg <- rbind(mut_plus_agg, mut_minus_agg)

#### create the plot

p.C.4 <- ggplot() + 
  geom_bar(data = mut_agg, aes(x = Group.2, y = log2(x)), 
           stat = "identity", position=position_dodge()) +
  facet_wrap(~direction) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  xlab("") + ylab("Plus/Minus \n(Log2)") +
  scale_fill_manual(values = repdomain_colors, guide = FALSE) 

# add plot format
p.C.4 <- p.C.4 + p_format + 
  theme(strip.background = element_blank(),
        strip.text.x = element_blank()) + ylim(-0.4, 0.4) 

p.C.4 # visualize


#### p.D.2 ####

# filter the data
d <- filter(df, phase != "async", replicate == "B", 
            product == "CPD", time_after_exposure == "12", phase == "late")

# create the plot 
p <- ggplot(d, aes(x = windows, y = log2(xr_ds))) + 
  geom_line(aes(color = sample_strand)) + 
  facet_grid(~repdomains) +
  xlab(windows_lab) + ylab(fr_xr_ds_lab) +
  scale_x_continuous(limits = c(-half_window, half_window), 
                     breaks = c(-half_window, 0, half_window), 
                     labels = c("-10", "0", "+10")) + 
  scale_color_manual(values = strand_colors) + 
  labs(color = "Strands")
p.D.2 <- p + p_format + ylim(-1.5, 1.5) + 
  geom_vline(xintercept = 0, color = "gray", linetype = "dashed") +
  theme(axis.title.x=element_blank(),
                                                axis.title.y=element_blank(),
                                                #axis.text.y=element_blank(),
                                                axis.text.x=element_text(hjust=c(0.3, 0.5, 0.7)))

p.D.2 # visualize



#### p.D.3 ####

mut_plus <- filter(d, windows > 0)
mut_plus_agg <- aggregate(x = mut_plus$xr_ds, by = 
                            list(mut_plus$dataset, mut_plus$repdomains, 
                                 mut_plus$sample_strand), FUN = "mean")
mut_plus_agg$direction <- "Right Replicating"
mut_minus <- filter(d, windows < 0)
mut_minus_agg <- aggregate(x = mut_minus$xr_ds, by = 
                             list(mut_minus$dataset, mut_minus$repdomains, 
                                  mut_minus$sample_strand), FUN = "mean")
mut_minus_agg$direction <- "Left Replicating"
mut_agg <- rbind(mut_plus_agg, mut_minus_agg)

#### create the plot 

p.D.3 <- ggplot() + 
  geom_bar(data = mut_agg, aes(x = Group.2, y = x, 
                               fill = Group.3), 
           stat = "identity", size = 1.5, 
           position=position_dodge2(padding = 0.05)) +
  facet_wrap(~direction) +
  xlab("Replication Domains") + 
  ylab("Repair \nRate") +
  scale_fill_manual(values = c("+" = "#0571b0", 
                               "-" = "#ca0020"), 
                    guide = FALSE) +
  labs(color = "Strands", fill = "") 

p.D.3 <- p.D.3 + p_format + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) + ylim(0, 2.5) # add plot format

#### p.D.4 ####

mut_casted <- dcast(d, dataset + repdomains + windows ~ sample_strand, 
                    value.var = "xr_ds")
mut_casted$plus_min <- mut_casted$"+" / mut_casted$"-" 
mut_plus <- filter(mut_casted, windows > 0)
mut_plus_agg <- aggregate(x = mut_plus$plus_min, by = 
                            list(mut_plus$dataset, mut_plus$repdomains), 
                          FUN = "mean")
mut_plus_agg$direction <- "Right Replicating"
mut_minus <- filter(mut_casted, windows < 0)
mut_minus_agg <- aggregate(x = mut_minus$plus_min, by = 
                             list(mut_minus$dataset, mut_minus$repdomains), 
                           FUN = "mean")
mut_minus_agg$direction <- "Left Replicating"
mut_agg <- rbind(mut_plus_agg, mut_minus_agg)

#### create the plot

p.D.4 <- ggplot() + 
  geom_bar(data = mut_agg, aes(x = Group.2, y = log2(x)), 
           stat = "identity", position=position_dodge()) +
  facet_wrap(~direction) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  xlab("") + ylab("Plus/Minus \n(Log2)") +
  scale_fill_manual(values = repdomain_colors, guide = FALSE) 

# add plot format
p.D.4 <- p.D.4 + p_format + 
  theme(strip.background = element_blank(),
        strip.text.x = element_blank()) + ylim(-0.4, 0.4) 

p.D.4 # visualize


#### p.E.2 ####

# filter the data
d <- filter(df, phase != "async", replicate == "B", 
            product == "CPD", time_after_exposure == "120", phase == "early")

# create the plot 
p <- ggplot(d, aes(x = windows, y = log2(xr_ds))) + 
  geom_line(aes(color = sample_strand)) + 
  facet_grid(~repdomains) +
  xlab(windows_lab) + ylab(fr_xr_ds_lab) +
  scale_x_continuous(limits = c(-half_window, half_window), 
                     breaks = c(-half_window, 0, half_window), 
                     labels = c("-10", "0", "+10")) + 
  scale_color_manual(values = strand_colors) + 
  labs(color = "Strands")
p.E.2 <- p + p_format + ylim(-1.5, 1.5) + 
  geom_vline(xintercept = 0, color = "gray", linetype = "dashed") +
  theme(axis.title.x=element_blank(),
                                                axis.title.y=element_blank(),
                                                #axis.text.y=element_blank(),
                                                axis.text.x=element_text(hjust=c(0.3, 0.5, 0.7)))

p.E.2 # visualize



#### p.E.3 ####

mut_plus <- filter(d, windows > 0)
mut_plus_agg <- aggregate(x = mut_plus$xr_ds, by = 
                            list(mut_plus$dataset, mut_plus$repdomains, 
                                 mut_plus$sample_strand), FUN = "mean")
mut_plus_agg$direction <- "Right Replicating"
mut_minus <- filter(d, windows < 0)
mut_minus_agg <- aggregate(x = mut_minus$xr_ds, by = 
                             list(mut_minus$dataset, mut_minus$repdomains, 
                                  mut_minus$sample_strand), FUN = "mean")
mut_minus_agg$direction <- "Left Replicating"
mut_agg <- rbind(mut_plus_agg, mut_minus_agg)

#### create the plot 

p.E.3 <- ggplot() + 
  geom_bar(data = mut_agg, aes(x = Group.2, y = x, 
                               fill = Group.3), 
           stat = "identity", size = 1.5, 
           position=position_dodge2(padding = 0.05)) +
  facet_wrap(~direction) +
  xlab("Replication Domains") + 
  ylab("Repair \nRate") +
  scale_fill_manual(values = c("+" = "#0571b0", 
                               "-" = "#ca0020"), 
                    guide = FALSE) +
  labs(color = "Strands", fill = "") 

p.E.3 <- p.E.3 + p_format + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) + ylim(0, 2.5) # add plot format

#### p.E.4 ####

mut_casted <- dcast(d, dataset + repdomains + windows ~ sample_strand, 
                    value.var = "xr_ds")
mut_casted$plus_min <- mut_casted$"+" / mut_casted$"-" 
mut_plus <- filter(mut_casted, windows > 0)
mut_plus_agg <- aggregate(x = mut_plus$plus_min, by = 
                            list(mut_plus$dataset, mut_plus$repdomains), 
                          FUN = "mean")
mut_plus_agg$direction <- "Right Replicating"
mut_minus <- filter(mut_casted, windows < 0)
mut_minus_agg <- aggregate(x = mut_minus$plus_min, by = 
                             list(mut_minus$dataset, mut_minus$repdomains), 
                           FUN = "mean")
mut_minus_agg$direction <- "Left Replicating"
mut_agg <- rbind(mut_plus_agg, mut_minus_agg)

#### create the plot

p.E.4 <- ggplot() + 
  geom_bar(data = mut_agg, aes(x = Group.2, y = log2(x)), 
           stat = "identity", position=position_dodge()) +
  facet_wrap(~direction) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  xlab("") + ylab("Plus/Minus \n(Log2)") +
  scale_fill_manual(values = repdomain_colors, guide = FALSE) 

# add plot format
p.E.4 <- p.E.4 + p_format + 
  theme(strip.background = element_blank(),
        strip.text.x = element_blank()) + ylim(-0.4, 0.4) 

p.E.4 # visualize


#### p.F.2 ####

# filter the data
d <- filter(df, phase != "async", replicate == "B", 
            product == "CPD", time_after_exposure == "120", phase == "late")

# create the plot 
p <- ggplot(d, aes(x = windows, y = log2(xr_ds))) + 
  geom_line(aes(color = sample_strand)) + 
  facet_grid(~repdomains) +
  xlab("Relative Position Centered \nat Initiation Zones (kb)") + 
  ylab(fr_xr_ds_lab) +
  scale_x_continuous(limits = c(-half_window, half_window), 
                     breaks = c(-half_window, 0, half_window), 
                     labels = c("-10", "0", "+10")) + 
  scale_color_manual(values = strand_colors) + 
  labs(color = "Strands")
p.F.2 <- p + p_format + ylim(-1.5, 1.5) + 
  geom_vline(xintercept = 0, color = "gray", linetype = "dashed") +
  theme(axis.title.y=element_blank(),
        #axis.text.y=element_blank(),
        axis.text.x=element_text(hjust=c(0.3, 0.5, 0.7)))

p.F.2 # visualize


#### p.F.3 ####

mut_plus <- filter(d, windows > 0)
mut_plus_agg <- aggregate(x = mut_plus$xr_ds, by = 
                            list(mut_plus$dataset, mut_plus$repdomains, 
                                 mut_plus$sample_strand), FUN = "mean")
mut_plus_agg$direction <- "Right Replicating"
mut_minus <- filter(d, windows < 0)
mut_minus_agg <- aggregate(x = mut_minus$xr_ds, by = 
                             list(mut_minus$dataset, mut_minus$repdomains, 
                                  mut_minus$sample_strand), FUN = "mean")
mut_minus_agg$direction <- "Left Replicating"
mut_agg <- rbind(mut_plus_agg, mut_minus_agg)

#### create the plot 

p.F.3 <- ggplot() + 
  geom_bar(data = mut_agg, aes(x = Group.2, y = x, 
                               fill = Group.3), 
           stat = "identity", size = 1.5, 
           position=position_dodge2(padding = 0.05)) +
  facet_wrap(~direction) +
  xlab("Replication Domains") + 
  ylab("Repair \nRate") +
  scale_fill_manual(values = c("+" = "#0571b0", 
                               "-" = "#ca0020"), 
                    guide = FALSE) +
  labs(color = "Strands", fill = "") 

p.F.3 <- p.F.3 + p_format + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) + ylim(0, 2.5) # add plot format

#### p.F.4 ####

mut_casted <- dcast(d, dataset + repdomains + windows ~ sample_strand, 
                    value.var = "xr_ds")
mut_casted$plus_min <- mut_casted$"+" / mut_casted$"-" 
mut_plus <- filter(mut_casted, windows > 0)
mut_plus_agg <- aggregate(x = mut_plus$plus_min, by = 
                            list(mut_plus$dataset, mut_plus$repdomains), 
                          FUN = "mean")
mut_plus_agg$direction <- "Right Replicating"
mut_minus <- filter(mut_casted, windows < 0)
mut_minus_agg <- aggregate(x = mut_minus$plus_min, by = 
                             list(mut_minus$dataset, mut_minus$repdomains), 
                           FUN = "mean")
mut_minus_agg$direction <- "Left Replicating"
mut_agg <- rbind(mut_plus_agg, mut_minus_agg)

#### create the plot

p.F.4 <- ggplot() + 
  geom_bar(data = mut_agg, aes(x = Group.2, y = log2(x)), 
           stat = "identity", position=position_dodge()) +
  facet_wrap(~direction) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  xlab("Replication Domains") + ylab("Plus/Minus \n(Log2)") +
  scale_fill_manual(values = repdomain_colors, guide = FALSE) 

# add plot format
p.F.4 <- p.F.4 + p_format + 
  theme(strip.background = element_blank(),
        strip.text.x = element_blank()) + ylim(-0.4, 0.4) 

p.F.4 # visualize




#### final plot ####

(p.A.1 + p.A.2 + (p.A.3 / p.A.4)) / (p.B.1 + p.B.2 + (p.B.3 / p.B.4)) / 
  (p.C.1 + p.C.2 + (p.C.3 / p.C.4)) / (p.D.1 + p.D.2 + (p.D.3 / p.D.4)) / 
  (p.E.1 + p.E.2 + (p.E.3 / p.E.4)) / (p.F.1 + p.F.2 + (p.F.3 / p.F.4)) + 
  plot_layout(guides = "collect") & 
  theme(legend.position = 'bottom') 

ggsave("~/Desktop/fig5_repB.png", width = 397, height = 410, units = "mm")
