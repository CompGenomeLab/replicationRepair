#### Packages and Libraries ####

library(stringr)
library(ggpubr)
library(ggplot2)
library(dplyr)
library(reshape2)
library(patchwork)
library(ggthemes)
library(grid)


#### Variables ####

# name of the sample csv file 
sample_csv <- paste("/home/azgarian/Documents/myprojects/replicationRepair/",
                    "final/final_reports_hg19_late_peaks_uv_repdomains_rep2",
                    "_windows_11_100.bed_ready.csv", 
                    sep = "")
sample_df <- read.csv( sample_csv )

type <- ""
edu_phase <- "late"
treatment <- "uv"

#### Default Plot Format ####

source("/home/azgarian/Documents/myprojects/replicationRepair/1_code/4_plot_format.R")


#### Fuctions ####

source("/home/azgarian/Documents/myprojects/replicationRepair/1_code/4_functions.R")


#### Main ####


for (myphase in c("late","early","async")) {

  
  if (edu_phase == "early"){ 
    
    edu <- "\nEarly EdU Peaks (bp)"
    edu_file <- "early_edu_rep2_cpd_" 
    
  } else if (myphase == "late"){ 
    
    edu <- "\nLate EdU Peaks (bp)"
    edu_file <- "late_edu_rep2_cpd_"   
    
  }
  
  if (treatment == "uv"){uv <- " UV treated"} else if (treatment == ""){uv <- ""}
  
  if (myphase == "early"){ 
    
    phase_name <- 'Early S Phase'
    myheight <- 18
    
  } else if (myphase == "late"){ 
    
    phase_name <- 'Late S Phase' 
    myheight <- 18
    
  } else if (myphase == "async"){ 
    
    myheight <- 9
  } 
  
  
  df_rr <- repair_rate( sample_df )
  df_rr_org <- window_numbering( df_rr, 4, 6 )
  df_rr_org <- domain_name( df_rr_org, 1 )
  df_rr_org$dataset <- "Peaks"
  
  df_rr_org$sample_strand <- factor(
    df_rr_org$sample_strand, levels = c("+","-"))
  
  # filtering for B.1
  pB1_data <- filter(df_rr_org, replicate == "A", 
                     product == "CPD", time_after_exposure == "12", 
                     phase == myphase)
  
  # for plot B.2
  pB2_data <- rr_boxplot( pB1_data ) 
  
  # for plot B.3
  pB3_data <- rr_boxplot_plus_minus( pB1_data ) 
  
  if (myphase == "async"){ } else {
  # filtering for C.1
  pC1_data <- filter(df_rr_org, replicate == "A", 
                     product == "CPD", time_after_exposure == "120", 
                     phase == myphase)
  
  # for plot C.2
  pC2_data <- rr_boxplot( pC1_data ) 
  
  # for plot C.3
  pC3_data <- rr_boxplot_plus_minus( pC1_data ) 
  }
  
  #### Plot A ####
  
  # plot A will be a drawing
  # we are creating an empty text for that part
  p.A <- wrap_elements(grid::textGrob(''))
  
  
  #### Plot B.1 ####
  
  # create the plot 
  p.B.1 <- ggplot(pB1_data, aes(x = windows, y = log2(xr_ds))) +
    geom_vline(xintercept = 0, color = "gray", linetype = "dashed") + 
    geom_line(aes(color = sample_strand)) + 
    facet_grid(~repdomains) +
    xlab(windows_lab) + ylab(fr_xr_ds_lab) +
    scale_x_continuous(limits = c(-6, 6), 
                       breaks = c(-6, 0, 6), 
                       labels = c("-500", "0", "+500")) + 
    scale_color_manual(values = strand_colors) + 
    #ylim(-1.5, 1.5) + 
    labs(color = "Strands")
  
  # adding and overriding the default plot format
  p.B.1 <- p.B.1 + p_format + 
    theme(axis.title.x=element_blank(),
          axis.text.x=element_text(hjust=c(0.3, 0.5, 0.7)))
  
  
  #### Plot B.2 ####
  
  # "Plus Minus" writing in the plot
  dat_text <- data.frame(
    label = c("Plus\nMinus", ""),
    direction = c("Left Replicating", "Right Replicating"),
    x = c(2, 2), y = c(0.8, 0.8))
  
  # create the plot 
  p.B.2 <- ggplot() + 
    geom_bar(data = pB2_data, aes(x = Group.2, y = x, 
                                  fill = Group.3), 
             stat = "identity", size = 1.5, 
             position=position_dodge2(padding = 0.05)) +
    facet_wrap(~direction) +
    xlab("Replication Domains") + 
    ylab("Repair\nRate (RR)") +
    scale_fill_manual(values = c("+" = "#0571b0", 
                                 "-" = "#ca0020"), 
                      guide = FALSE) +
    labs(color = "Strands", fill = "") +
    #scale_y_continuous(breaks = c(0, 1, 2),
    #                   limits = c(0, 2)) +
    geom_text( data = dat_text, mapping = aes(x = x, y = y, label = label), 
               angle = 90, colour = "white", size = 2.8 ) 
  
  # adding and overriding the default plot format
  p.B.2 <- p.B.2 + p_format + 
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x = element_blank()) 
  
  
  #### Plot B.3 ####
  
  # create the plot 
  p.B.3 <- ggplot() + 
    geom_bar(data = pB3_data, aes(x = Group.2, y = x), 
             stat = "identity", position=position_dodge()) +
    facet_wrap(~direction) +
    geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
    xlab("") + ylab(expression(RR[p] - RR[m])) +
    #scale_y_continuous(breaks = c(-.2, 0, .2),
    #                   limits = c(-.25, .25)) +
    scale_fill_manual(values = repdomain_colors, guide = FALSE) 
  
  # adding and overriding the default plot format
  p.B.3 <- p.B.3 + p_format + 
    theme(strip.background = element_blank(),
          strip.text.x = element_blank())   
  
  if (myphase == "async"){ 
    
    p.B.2.3 <- (p.B.2 / p.B.3)   
    
    layout <- "
    BBBCCCD
    "
    p_comb <- p.B.1 + p.B.2.3 + grid::textGrob('CPD\n12 min.\nAsynchronized', 
                                               rot = -90, gp=gpar(fontsize=12), 
                                               y = unit(.55, "npc")) 
    
  } else {
    
  
  #### Plot C.1 ####
  
  # create the plot 
  p.C.1 <- ggplot(pC1_data, aes(x = windows, y = log2(xr_ds))) + 
    geom_vline(xintercept = 0, color = "gray", linetype = "dashed") +
    geom_line(aes(color = sample_strand)) + 
    facet_grid(~repdomains) +
    xlab(paste("Position Relative to", uv, edu, sep="")) + 
    ylab(fr_xr_ds_lab) +
    scale_x_continuous(limits = c(-6, 6), 
                       breaks = c(-6, 0, 6), 
                       labels = c("-500", "0", "+500")) + 
    scale_color_manual(values = strand_colors) + 
    #ylim(-1.5, 1.5) + 
    labs(color = "Strands")
  
  # adding and overriding the default plot format
  p.C.1 <- p.C.1 + p_format + 
    theme(axis.text.x=element_text(hjust=c(0.3, 0.5, 0.7)))
  
  
  #### Plot C.2 ####
  
  # create the plot 
  p.C.2 <- ggplot() + 
    geom_bar(data = pC2_data, aes(x = Group.2, y = x, 
                                  fill = Group.3), 
             stat = "identity", size = 1.5, 
             position=position_dodge2(padding = 0.05)) +
    facet_wrap(~direction) +
    xlab("Replication Domains") + 
    ylab("Repair\nRate (RR)") +
    scale_fill_manual(values = c("+" = "#0571b0", 
                                 "-" = "#ca0020"), 
                      guide = FALSE) +
    #scale_y_continuous(breaks = c(0, 1, 2),
    #                   limits = c(0, 2)) +
    labs(color = "Strands", fill = "") 
  
  # adding and overriding the default plot format
  p.C.2 <- p.C.2 + p_format + 
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x = element_blank()) 
  
  
  #### Plot C.3 ####
  
  # create the plot 
  p.C.3 <- ggplot() + 
    geom_bar(data = pC3_data, aes(x = Group.2, y = x), 
             stat = "identity", position=position_dodge()) +
    facet_wrap(~direction) +
    geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
    xlab("Replication Domains") + 
    ylab(expression(RR[p] - RR[m])) +
    scale_fill_manual(values = repdomain_colors, guide = FALSE) 
    #scale_y_continuous(breaks = c(-.2, 0, .2),
    #                   limits = c(-.25, .25))
  
  # adding and overriding the default plot format
  p.C.3 <- p.C.3 + p_format + 
    theme(strip.background = element_blank(),
          strip.text.x = element_blank()) 
  
  
  #### Combining Plots with Patchwork ####
  
  p.C.2.3 <- (p.C.2 / p.C.3)
  p.B.2.3 <- (p.B.2 / p.B.3) 
  
  layout <- "
    BBBCCCD
    EEEFFFG
    "
  
  p_comb <-  p.B.1 + p.B.2.3 + grid::textGrob(paste('CPD\n12 min.\n', phase_name, sep=""), 
                                              rot = -90, gp=gpar(fontsize=12), 
                                              y = unit(.55, "npc")) + 
    p.C.1 + p.C.2.3 + grid::textGrob(paste('CPD\n120 min.\n', phase_name, sep=""), 
                                     rot = -90, gp=gpar(fontsize=12), 
                                     y = unit(.62, "npc")) 
  

  
  }
  
  p_comb + plot_layout(design = layout, guides = "collect") & 
    theme(plot.tag = element_text(size = 12, face="bold"),
          legend.position = 'bottom', 
          plot.title = element_text(hjust = -0.2, vjust = 5, 
                                    size = 12, face="bold"))
  
  ggsave(paste("~/Desktop/repairRevision2/edu_figs/", edu_file, treatment, "_", 
               myphase, type, ".png", sep=""), 
         width = 22, height = myheight, units = "cm")

}
  