#### library ####

library(ggplot2)
library(stringr)
library(dplyr)

#### prepare data ####

sourcePath <- paste("~/Documents/myprojects/replicationRepair/", 
                    "1_code/R/", sep = "") 

date <- "[2020.02.18.23:19]"

fr_name <- "hi_rfd_repdomains_mutations_combined"

#### set variables ####

dateout <- Sys.Date()
dateout <- format(dateout, format = "[%Y.%m.%d]")

dataInfoPath <- paste("~/Documents/myprojects/replicationRepair/", 
                 "4_output/gitignore/1_TextforPlotting/", sep = "")

#### import and rearrange data ####

mut <- read.table(paste(dataInfoPath, date, fr_name, ".bed", sep = ""))

names(mut) <- c("chr", "start_pos", "end_pos", "name", "score", "strands", "count")

Windows = data.frame(str_split_fixed(mut$name, "_", -1))

mut$window_number = as.numeric(levels(Windows[ , ncol(Windows)]
))[Windows[ , ncol(Windows)]] -101

mut$repdomains <- Windows[,2]

mut$name <- "hi_rfd"

mut$repdomains <- factor(mut$repdomains, levels = c("UTZ", "ERD", "DTZ", "LRD"))

#### plot ####

#### filter the data ####

#d <- filter(mut, repdomains == "DTZ") # filter

#### add plot format #### 

source(paste(sourcePath, "4_plot_format.R", sep = ""))

#### create the plot ####

p <- ggplot(mut) + 
  geom_line(aes(x = window_number, y = count, 
                color = strands)) +  
  facet_grid(~score~repdomains, labeller = 
               labeller(score = dataset_strand_labs)) +
  xlab(windows_lab) + ylab("Mutation Count per Region") +
  scale_x_continuous(limits = c(-100, 100), 
                     breaks = c(-90, 0, 90), 
                     labels = c("-1 Mb", "high RFDs", "+1 Mb")) + 
  scale_color_manual(values = strand_colors) + 
  labs(color = "Strands") +
  labs(title = paste("Mutation Profile of high RFDs"))

p <- p + p_format # add plot format

p # visualize

#### save the plot ####

figurePath <- paste("~/Documents/myprojects/replicationRepair/", 
                    "4_output/gitignore/Target_mutations/", sep = "")

figureName <- paste(dateout, 
                    "Mutation_Profile_of_hi_rfd_with_repdomains.pdf", 
                    sep = "")

ggsave(path = figurePath, filename = figureName, 
       width = 297, height = 210, units = "mm")

figurePNG <- sub(".pdf", ".png", figureName)

ggsave(path = figurePath, filename = figurePNG, 
       width = 297, height = 210, units = "mm")

#### save the figure info ####

dataInfo <- paste(date, "final_report_", fr_name, "_info.TXT", sep = "")

figureInfo <- "3_target_mutations_intersected_Inzones.R"

source(paste(sourcePath, "4_figure_info.R", sep = ""))

rm(mut, p, Windows, dateout)
