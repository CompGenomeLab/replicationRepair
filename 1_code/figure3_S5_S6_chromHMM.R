#### Packages and Libraries ####

library(ggpubr)
library(ggplot2)
library(stringr)
library(reshape2)
library(dplyr)
library(ggthemes)
library(patchwork)


#### Variables ####

rep <- "A" # replicate
dprod <- "CPD" # damage product
if (rep == "A" && dprod == "CPD"){ fig_name = "~/Desktop/fig3.png" 
} else if (rep == "B" && dprod == "CPD"){ fig_name = "~/Desktop/supfig5.png"
} else if (rep == "A" && dprod == "64_PP"){ fig_name = "~/Desktop/supfig6.png" 
} 

# name of the sample csv file
sample_csv <- paste("~/Documents/myprojects/replicationRepair/3_output/",
                    "gitignore/1_TextforPlotting/[2020.04.07]final_report_", 
                    "chromhmm_windows_chr_ready.csv", 
                    sep = "")


#### Default Plot Format ####

source("4_plot_format.R")


#### Fuctions ####

source("4_functions.R")


#### Main ####

# for both plot A and plot B
sample_df <- read.csv( sample_csv )
df_rr <- repair_rate( sample_df )

# for plot A
pA_df_rr <- chrState_naming( df_rr, chr_states, general_states, 
                             chrState2generalState )

# for plot B
df_rr_ear_la <- rrEarly_rrLate( df_rr )
pB_df_rr_ear_la <- chrState_naming( df_rr_ear_la, chr_states, general_states, 
                                    chrState2generalState )


#### Filtering Samples ####

# for plot A
pA_filt <- filter(pA_df_rr, xr_ds != "NaN", xr_ds != "Inf", 
                  dataset_strand != "DTZ", dataset_strand != "UTZ", 
                  product == dprod, replicate == rep, phase != "async", 
                  time_after_exposure == "12")
pA_data <- aggregate(x = pA_filt$xr_ds, by = list(pA_filt$dataset, 
                                                  pA_filt$dataset_strand, 
                                                  pA_filt$phase, 
                                                  pA_filt$states, 
                                                  pA_filt$chromosomes), 
                     FUN = "mean")


# for plot B
pB_filt <- filter(pB_df_rr_ear_la, ear_la != "NaN", ear_la != "Inf", 
                  dataset_strand != "DTZ", dataset_strand != "UTZ", 
                  product == dprod, replicate == rep,
                  time_after_exposure == "12")
pB_data <- aggregate(x = pB_filt$ear_la, by = list(pB_filt$dataset, 
                                                   pB_filt$dataset_strand, 
                                                   pB_filt$states, 
                                                   pB_filt$chromosomes), 
                     FUN = "mean")


#### Plot A ####

# create the plot
p.A <- ggplot(pA_data, aes(x = Group.1, y = log2(x))) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_boxplot(aes(fill = factor(Group.4)), outlier.shape = NA, lwd=0.1) +
  facet_grid(~Group.3~Group.2, labeller = labeller(Group.3 = phase_labs)) +
  xlab("") + ylab("Repair Rate (RR) (log2)") +
  scale_fill_manual(name = "Chromatin States", values = state_colors) +
  ylim(-2,3) +
  guides(size = FALSE) 

# adding and overriding the default plot format
p.A <- p.A + p_format + 
  theme(panel.border=element_rect(size=1, fill = NA),
        legend.position = "right", 
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())


#### Plot B ####

# create the plot
p.B <- ggplot(pB_data, aes(x = Group.1, y = log2(x))) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_boxplot(aes(fill = factor(Group.3)), outlier.shape = NA, lwd=0.1) +
  facet_wrap(~Group.2) +
  xlab(chrState_lab) + ylab(expression(RR[E] / RR[L] (log2))) +
  scale_fill_manual(name = "Chromatin States", values = state_colors) +
  ylim(-1,1) +
  guides(size = FALSE) 

# adding and overriding the default plot format
p.B <- p.B + p_format +
  theme(panel.border=element_rect(size=1, fill = NA),
        strip.text = element_blank(),
        legend.position = "right", 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.95))


#### Combining Plots with Patchwork ####

p.A / p.B + 
  plot_layout(guides = "collect") & 
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size = 12, face="bold"),
        legend.margin=margin())

ggsave(fig_name, width = 22, height = 18, units = "cm")
