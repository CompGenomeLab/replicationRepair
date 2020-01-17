#### library ####

library(ggplot2)
library(stringr)

#### set variables ####

dateout <- Sys.Date()
dateout <- format(dateout, format = "[%Y.%m.%d]")

mut_dir <- paste("~/Documents/My_Projects/Project_Repair_Replication/", 
                 "Results/1_TextforPlotting/", sep = "")

#### import and rearrange data ####

mut_raw <- read.table(paste(mut_dir, 
                            "[2019.09.25]Inzones_mutations_combined.bed", 
                            sep = ""))

mut_org <- mut_raw[ , c(2,3,4,5)]

names(mut_org) <- c("name", "counts", "window.counts", "strands")

Windows = data.frame(str_split_fixed(mut_org$name, "_", -1))

mut_org$window_number = as.numeric(levels(Windows[ , ncol(Windows)]
                                          ))[Windows[ , ncol(Windows)]] -101

mut_org$name <- "Initiation_Zones"

mut_reorg = mut_org[ , c("name", "window_number", "strands", "counts")]

#### plot ####

p <- ggplot(mut_reorg) + 
  geom_line(aes(x = window_number, y = counts, color = strands)) + 
  geom_smooth(method = "loess", aes(x = window_number, y = counts, 
                                    color = strands), se = FALSE) + 
  xlab("Windows") + ylab("Counts") +
  scale_color_manual(name = "Strands", values = c("blue", "red")) +
  scale_x_continuous(limits = c(-100, 100), 
                     breaks = c(-100, 0, 100), 
                     labels = c("-10 kb", "0", "+10 kb")) + 
  labs(title = paste("Mutation Profiles of Initiation Zones"), 
       subtitle = paste("20 kb sliding windows at 100 bp intervals", sep = ""))

#### theme ####

p <- p + theme_light() +
          theme(axis.title.x = element_text(size = 20, face = "bold"),
                axis.title.y = element_text(size = 20, face = "bold"),
                axis.text.x = element_text(size = 18, vjust = 0.6),
                axis.text.y = element_text(size = 18, vjust = 0.1),
                strip.text.x = element_text(size = 18),
                strip.text.y = element_text(size = 18, angle = 360),
                legend.title = element_text(size = 18, face = "bold"),
                legend.text = element_text(size = 18))
          
setwd(paste("~/Documents/My_Projects/Project_Repair_Replication/", 
            "Results/Target_mutations/", sep = ""))

ggsave(filename = paste(dateout, "Mutation_Profile_of_InitiationZones_200", 
                        "_sliding_windows_at_100bp_intervals.pdf", sep = ""), 
       width = 297, height = 210, units = "mm")

rm(mut_org, mut_raw, mut_reorg, p, Windows, dateout, mut_dir)
