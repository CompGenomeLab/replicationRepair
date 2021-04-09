#### Packages and Libraries ####

library(stringr)
library(ggpubr)
library(ggplot2)
library(dplyr)
library(reshape2)
library(patchwork)
library(ggthemes)


#### Variables ####

rep <- "B" 
if (rep == "A"){ 
  fig_name = "~/Desktop/fig2.svg" 
  ann_text_ERD <- data.frame(phase = 0.7, xr_ds = 4.8,
                             dataset = factor("ERD",levels = c("ERD","LRD")),
                             product = factor("64_PP", 
                                              levels = c("64_PP", "CPD")))
  lab_erd <- "p ="
  ann_text_LRD <- data.frame(phase = 0.9, xr_ds = 4.8,
                             dataset = factor("LRD",levels = c("ERD","LRD")),
                             product = factor("64_PP", 
                                              levels = c("64_PP", "CPD")))
  lab_lrd <- "p ="
} else if (rep == "B"){ 
  fig_name = "~/Desktop/supfig4.svg" 
  ann_text_ERD <- data.frame(phase = 0.7, xr_ds = 4.8,
                             dataset = factor("ERD",levels = c("ERD","LRD")),
                             product = factor("64_PP", 
                                              levels = c("64_PP", "CPD")))
  lab_erd <- ""
  ann_text_LRD <- data.frame(phase = 0.7, xr_ds = 4.7,
                             dataset = factor("LRD",levels = c("ERD","LRD")),
                             product = factor("64_PP", 
                                              levels = c("64_PP", "CPD")))
  lab_lrd <- "p ="
} 

# name of the sample csv file for plot A
sample_csv_pA <- paste("~/Documents/myprojects/replicationRepair/3_output/",
                        "gitignore/1_TextforPlotting/", 
                        "[2019.10.31]final_report_ERD_LRD_windows_ready.csv", 
                        sep = "")

# name of the sample csv file for plot B
sample_csv_pB <- paste("~/Documents/myprojects/replicationRepair/3_output/",
                       "gitignore/1_TextforPlotting/", 
                       "[2019.10.31]final_report_repdomain_ready.csv", 
                       sep = "")

# path of the default plot format and functions
sourcePath <- "~/Documents/myprojects/replicationRepair/1_code/r/"


#### Default Plot Format ####

source(paste(sourcePath, "4_plot_format.R", sep = ""))


#### Fuctions ####

source(paste(sourcePath, "4_functions.R", sep = ""))


#### Main ####

# for plot A
pA_sample_df <- read.csv( sample_csv_pA )
pA_df_rr <- repair_rate( pA_sample_df )
pA_df_rr_org <- window_numbering( pA_df_rr, 4, 101 )
pA_df_rr_org$dataset <- gsub("_.*", "", pA_df_rr_org$dataset)

# for plot B
pB_sample_df <- read.csv( sample_csv_pB )
pB_sample_df <- rmv_low_counted_regions( pB_sample_df )
pB_df_rr <- repair_rate( pB_sample_df )

#### Filtering Samples ####

pA1_data <- filter(pA_df_rr_org, phase != "async", dataset == "ERD", 
                   replicate == rep, time_after_exposure == "12")

pA2_data <- filter(pA_df_rr_org, phase != "async", dataset == "LRD", 
                   replicate == rep, time_after_exposure == "12")

pB_data = filter(pB_df_rr, dataset == "ERD" | dataset == "LRD", 
                 replicate == rep, phase != "async", 
                 time_after_exposure == "12")


#### Plot A.1 ####

# create the plot
p.A.1 <- ggplot(pA1_data, aes(x = windows, y = log2(xr_ds))) + 
  geom_hline(yintercept = 0, linetype="dashed", color="red") +
  geom_line(aes(color = phase, linetype = sample_strand)) +
  facet_grid(~product~time_after_exposure~dataset,
             labeller = labeller(product = product_labs, 
                                 time_after_exposure = taex_labs)) +
  xlab("Relative Position (kb)") + ylab(fr_xr_ds_lab) +
  scale_x_continuous(limits = c(-101, 101), 
                     breaks = c(-101, 0, 101), 
                     labels = c("-1000", "0", "+1000")) + 
  scale_color_manual(name = "Phase", 
                     label = c("Early Phase", "Late Phase", "Asyncronized"), 
                     values = phase_colors) + 
  labs(color = "Phase", linetype = "Strands") +
  ylim(-1.6, 1) +
  guides(linetype = FALSE)

# adding and overriding the default plot format
p.A.1 <- p.A.1 + p_format + 
  theme(strip.text.y = element_text(size=0, margin = margin(0, 0, 0, 0)),
        legend.position = "right",
        axis.text.x = element_text(hjust=c(0.1, 0.5, 0.9)))


#### Plot A.2 ####

# create the plot
p.A.2 <- ggplot(pA2_data, aes(x = windows, y = log2(xr_ds))) +
  geom_hline(yintercept = 0, linetype="dashed", color="red") +
  geom_line(aes(color = phase, linetype = sample_strand)) +
  facet_grid(~product~time_after_exposure~dataset,
             labeller = labeller(product = product_labs, 
                                 time_after_exposure = taex_labs)) +
  xlab("Relative Position (kb)") + ylab(fr_xr_ds_lab) +
  scale_x_continuous(limits = c(-101, 101), 
                     breaks = c(-101, 0, 101), 
                     labels = c("-1000", "0", "+1000")) + 
  scale_color_manual(name = "Phase", 
                     label = c("Early Phase", "Late Phase", "Asyncronized"), 
                     values = phase_colors) + 
  labs(color = "Phase", linetype = "Strands") +
  ylim(-1.6, 1) +
  guides(linetype = FALSE)

# adding and overriding the default plot format
p.A.2 <- p.A.2 + p_format + 
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        strip.text.y = element_text(size=0, margin = margin(0, 0, 0, 0)),
        legend.position = "right",
        axis.text.x = element_text(hjust=c(0.1, 0.5, 0.9)))


#### Plot B ####

# create the plot
p.B <- ggplot(pB_data, aes(x = phase, y = log2(xr_ds))) + 
  geom_boxplot(aes(fill = phase), outlier.shape = NA) +
  geom_hline(yintercept = 0, linetype="dashed", color="red") +
  facet_grid(product~time_after_exposure~dataset, 
             labeller = labeller(product = product_labs, 
                                 time_after_exposure = taex_labs, 
                                 replicate = rep_labs)) + 
  xlab(phase_lab) + ylab("") +
  scale_x_discrete(labels = c("Early S", "Late S")) +
  scale_fill_manual(name="", 
                    values = phase_colors) +
  guides(fill = FALSE) 

# adding and overriding the default plot format
# p values calculated
p.B <- p.B + p_format + 
  theme(axis.title.y=element_blank()) +
  stat_compare_means( comparisons = list( c("early", "late") ), 
                      label.y = 2) +
  geom_text(data = ann_text_ERD, label = lab_erd) +
  geom_text(data = ann_text_LRD, label = lab_lrd)


#### Combining Plots with Patchwork ####
p.A.1 + p.A.2 + p.B +
  plot_layout(guides = "collect") &
  plot_annotation(tag_levels = 'fig2') &
  theme(plot.tag = element_text(size = 12, face="bold"),
        legend.position = 'bottom', 
        plot.title = element_text(hjust = -0.2, vjust = 5, 
                                  size = 12, face="bold"))

ggsave(fig_name, width = 22, height = 18, units = "cm")



