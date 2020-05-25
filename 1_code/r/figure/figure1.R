library(stringr)
library(ggpubr)
library(ggplot2)
library(dplyr)
library(patchwork)
library(ggthemes)

#### erd early late plot ####

#### prepare data

sourcePath <- paste("~/Documents/myprojects/replicationRepair/1_code/R/", 
                    sep = "") 

date <- "[2019.10.31]"

fr_name <- "ERD_LRD_windows"

source(paste(sourcePath, "2_report_sub_dfs.R", sep = ""))

#### set df, variables and rearrange 

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

#### plot 

#### filter the data 

d <- filter(df, phase == "early" | phase == "late", dataset == "ERD", 
            replicate == "A")

#### add plot format

source(paste(sourcePath, "4_plot_format.R", sep = ""))

#### create the plot 

p.A.1.1 <- ggplot(d, aes(x = windows, y = log2(xr_ds))) + 
  geom_line(aes(linetype = sample_strand, color = phase)) +
  geom_hline(yintercept = 0, linetype="dashed", color="red") +
  facet_grid(~product~time_after_exposure~dataset,
             labeller = labeller(product = product_labs, 
                                 time_after_exposure = taex_labs)) +
  xlab("Relative Position (kb)") + ylab(fr_xr_ds_lab) +
  scale_x_continuous(limits = c(-half_window, half_window), 
                     breaks = c(-half_window, 0, half_window), 
                     labels = c("-1000", "0", "+1000")) + 
  scale_color_manual(name = "Phase", 
                     label = c("Early Phase", "Late Phase", "Asyncronized"), 
                     values = phase_colors) + 
  labs(color = "Phase", linetype = "Strands") +
  guides(linetype = FALSE)

p.A.1.1 <- p.A.1.1 + p_format + ylim(-1.6, 1) +
  theme(strip.text.y = element_blank(),
        legend.position = "right")

p.A.1.1 # visualize

#### filter the data

d <- filter(df, phase == "early" | phase == "late", dataset == "LRD", 
            replicate == "A")

#### add plot format

source(paste(sourcePath, "4_plot_format.R", sep = ""))

#### create the plot

p.A.1.2 <- ggplot(d, aes(x = windows, y = log2(xr_ds))) + 
  geom_line(aes(linetype = sample_strand, color = phase)) +
  geom_hline(yintercept = 0, linetype="dashed", color="red") +
  facet_grid(~product~time_after_exposure~dataset,
             labeller = labeller(product = product_labs, 
                                 time_after_exposure = taex_labs)) +
  xlab("Relative Position (kb)") + ylab(fr_xr_ds_lab) +
  scale_x_continuous(limits = c(-half_window, half_window), 
                     breaks = c(-half_window+10, 0, half_window-10), 
                     labels = c("-1000", "0", "+1000")) + 
  scale_color_manual(name = "Phase", 
                     label = c("Early Phase", "Late Phase", "Asyncronized"), 
                     values = phase_colors) + 
  labs(color = "Phase", linetype = "Strands") +
  guides(linetype = FALSE)

p.A.1.2 <- p.A.1.2 + p_format + ylim(-1.6, 1) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        strip.text.y = element_blank(),
        legend.position = "right")

p.A.1.2 # visualize

#### samples diff ####

#### prepare data 

sourcePath <- paste("~/Documents/myprojects/replicationRepair/1_code/R/", 
                    sep = "") 

date <- "[2019.10.31]"

fr_name <- "repdomain"

source(paste(sourcePath, "2_report_sub_dfs.R", sep = ""))

#### set df, variables and rearrange 

dateout <- Sys.Date()
dateout <- format(dateout, format = "[%Y.%m.%d]")
options(scipen = 999)

df <- fr_xr_ds

#### plot

#### filter the data 

d = filter(df, dataset == "ERD" | dataset == "LRD", replicate == "A", 
           phase != "async")

#### add plot format 

source(paste(sourcePath, "4_plot_format.R", sep = ""))

phase_labs <- c("Asyncronized", "Early Phase", "Late Phase")
names(phase_labs) <- c("async", "early", "late")

#### create the plot 

p.A.1.3 <- ggplot(d, aes(x = phase, y = log2(xr_ds))) + 
  geom_boxplot(aes(fill = phase), outlier.shape = NA) +
  geom_hline(yintercept = 0, linetype="dashed", color="red") +
  facet_grid(product~time_after_exposure~dataset, 
             labeller = labeller(product = product_labs, 
                                 time_after_exposure = taex_labs, 
                                 replicate = rep_labs)) + 
  xlab(phase_lab) + ylab("") +
  scale_x_discrete(labels = c("Early", "Late")) +
  scale_fill_manual(name="", 
                    values = phase_colors) +
  guides(fill = FALSE) 

p.A.1.3 <- p.A.1.3 +  p_format +  ylim(-2.5, 2.7) +
  stat_compare_means(
  label.sep = "\n",
  label.y = 2.3,
  label.x.npc = c("left")) +
  #stat_compare_means(method = "t.test", label.y = 1.8, label.x = 0.6) +
  theme(axis.title.y=element_blank()
        #strip.text.y = element_blank(),
        #strip.text.x = element_blank()
        #axis.text.y=element_blank(),
        #axis.ticks.y=element_blank()
        ) 

p.A.1.3  # visualize

((p.A.1.1 + p.A.1.2 + p.A.1.3) + plot_layout(guides = "collect") & 
    theme(legend.position = 'bottom') ) 

ggsave("~/Desktop/fig5_sub2.png", width = 380, height = 210, units = "mm")
