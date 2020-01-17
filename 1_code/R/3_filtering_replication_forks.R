#### library ####

library(ggplot2)
library(dplyr)
library(stringr)
library(reshape2)
library(ggnewscale)

#### set variables ####

options(scipen = 999)
dateout <- Sys.Date()
dateout <- format(dateout, format = "[%Y.%m.%d]")

setwd(paste("~/Documents/My_Projects/Project_Repair_Replication/",
            "Results/1_TextforPlotting/", sep = ""))

#### import and rearrange data ####

genome <- read.table("[2019.12.09]final_report_20kb_windows.txt", 
                     header = FALSE)

names(genome) <- c("chromosomes", "start_position", "end_position", 
                   "counts", "sample_names", 
                   "file_names", "layout", "cell_line", "product", "method", 
                   "uv_exposure", "treatment", "phase", "time_after_exposure", 
                   "replicate", "project", "sample_source", "sample_strand", 
                   "mapped_reads", "RPKM")

genome_early2late = filter(genome, 
                           sample_names == "Hela_NA_DNAseq_NA_early_NA_NA" | 
                             sample_names == "Hela_NA_DNAseq_NA_late_NA_NA", 
                           phase != "async" ) 

genome_early2late <- dcast(genome_early2late, chromosomes + start_position + 
                     end_position + sample_strand ~ phase, value.var = "RPKM")

genome_early2late$ear_la <- genome_early2late$early / genome_early2late$late

#genome_early2late = select(genome_early2late, -c(early, late))


inZones <- read.table(paste("[2019.12.09]InitiationZones", 
                            "_genome_hg19_20kb_windows.txt", sep = ""),
                      header = FALSE)

inZones$position <- NA

inZones <- within(inZones, position[V7 == "."] <- 3)

inZones <- within(inZones, position[V7 == "Initiation_Zones"] <- 1)

inZones = inZones[ ,c("V1","V2","V3","position")]

names(inZones) <- c("chromosomes", "start_position", "end_position", 
                    "InZones")

genome_final <- merge(genome_early2late, inZones, 
                     by.x=c("chromosomes", "start_position", "end_position"), 
                     by.y=c("chromosomes", "start_position", "end_position"), 
                     all.x=TRUE)

#### plot ####

d <- filter(genome, chromosomes == "chr1") #early > .2, ear_la > 1) 

chr_chosen <- d[1,1]
chr_name <- sub("chr", "Chromosome ", chr_chosen)

p <- ggplot(d) + 
  geom_point(aes(x = start_position, y = early, color = "Early Phase")) + 
  geom_point(aes(x = start_position, y = late, color = "Late Phase")) + 
  scale_color_manual(name = "Input", values = c("red", "blue")) +
  new_scale("color") +
  xlab(paste(chr_name, " (Mbp)", sep = "")) + ylab("RPKM") +
  ylim(0, 1.5) +
  scale_x_continuous(limits = c(0, 30000000), 
                     breaks = c(0, 5000000, 10000000, 15000000, 20000000, 
                                25000000, 30000000), 
                     labels = c("0", "5", "10", "15", "20", "25", "30")) + 
  guides(size = FALSE, alpha = FALSE) +
  geom_point(aes(x = start_position, y = ear_la, color = ear_la)) +
  geom_point(aes(x = start_position, y = InZones)) +
  scale_colour_gradient(name = "Phase", 
                        low = "blue", high = "red", breaks=c(0.75, 1.25), 
                        labels=c("Late", "Early"),
                        limits=c(0.75, 1.25)) 

#### theme ####

p <- p + theme_light() +
  theme(axis.title.x = element_text(size = 20, face = "bold"),
        axis.title.y = element_text(size = 20, face = "bold"),
        axis.text.x = element_text(size = 18, vjust = 0.6),
        axis.text.y = element_text(size = 18, vjust = 0.1),
        strip.text.x = element_text(size = 18),
        strip.text.y = element_text(size = 18, angle = 360),
        legend.title = element_text(size = 18, vjust = 1, face = "bold"),
        legend.text = element_text(size = 14),
        legend.position = "bottom",
        legend.box = "vertical")

p

setwd("~/Desktop")

ggsave(filename = "region.pdf", width = 297, height = 210, units = "mm")

#### export the positions as a bed file ####

setwd("~/Desktop")

genomeBed <- filter(genome_final, early > .2, ear_la > 1.2) 

genomeBed$name <- "Fork_Positions"

genomeBed$score <- "0"

genomeBed$strand <- "."

genomeBed <- genomeBed[ , c(1,2,3,9,10,11)]

write.table(genomeBed, file = "replication_fork.bed", sep = "\t", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)

d <- filter(genome_final, chromosomes == "chr1", phase == "early", InZones == "0.5") 

ggplot(d, aes(x = RPKM)) + 
  geom_histogram(aes(fill = phase)) + 
  facet_grid(product~time_after_exposure~replicate~InZones) +
  xlim(0,0.7)

setwd("~/Desktop")

ggsave(filename = "notinzones.pdf", width = 297, height = 210, units = "mm")

