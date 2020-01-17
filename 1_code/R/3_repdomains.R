#### library ####

library(ggplot2)
library(stringr)
library(ggpubr)

#### prepare data ####

sourcePath <- paste("~/Documents/My_Projects/Project_Repair_Replication/", 
                    "Scripts/R/", sep = "") 

date <- "[2019.10.31]"

fr_name <- "repdomain"

source(paste(sourcePath, "2_report_sub_dfs.R", sep = ""))

#### set df, variables and rearrange ####

dateout <- Sys.Date()
dateout <- format(dateout, format = "[%Y.%m.%d]")
options(scipen = 999)

df <- fr_xr_ds

df$chromosomes <- factor(df$chromosomes, levels = c("chr1", "chr2", "chr3", 
                                                  "chr4", "chr5", "chr6", 
                                                  "chr7", "chr8", "chr9", 
                                                  "chr10", "chr11", "chr12", 
                                                  "chr13", "chr14", "chr15", 
                                                  "chr16", "chr17", "chr18", 
                                                  "chr19", "chr20", "chr21", 
                                                  "chr22", "chrX"))

#### plot ####

#### filter the data ####

d = filter(df, dataset == "ERD" | dataset == "LRD", replicate == "A", 
           product == "CPD", time_after_exposure == "12")

#### add plot format #### 

source(paste(sourcePath, "4_plot_format.R", sep = ""))

#### create the plot ####

p <- ggplot(d, aes(x = chromosomes, y = log2(xr_ds))) + 
  geom_boxplot(aes(fill = dataset), outlier.shape = NA) +
  geom_hline(yintercept = -1, linetype="dashed", color = "darkred") +
  geom_hline(yintercept = 1, linetype="dashed", color = "darkgreen") +
  facet_grid(phase~time_after_exposure~replicate~product,
             labeller = labeller(product = product_labs, 
                                 time_after_exposure = taex_labs, 
                                 replicate = rep_labs, phase = phase_labs)) +
  xlab("Chromosomes") + ylab(fr_xr_ds_lab) +
  scale_fill_manual(name="", 
                    values = repdomain_colors) +
  ylim(-3,2)


p <- p +  p_format + theme(axis.text.x = 
                             element_text(size = 12, vjust = 0.6, angle = 60))

p # visualize

#### save the plot ####

figurePath <- paste("~/Documents/My_Projects/Project_Repair_Replication/", 
                    "Results/Repdomains/", sep = "")

figureName <- paste(dateout, "Repair_over_Damage_Value_of_", 
                    fr_name, "_for_Every_Chromosome_CPD_12.pdf", sep = "")

ggsave(path = figurePath, filename = figureName, 
       width = 297, height = 210, units = "mm")

figurePNG <- sub(".pdf", ".png", figureName)

ggsave(path = figurePath, filename = figurePNG, 
       width = 297, height = 210, units = "mm")

#### save the figure info ####

dataInfo <- paste(date, "final_report_", fr_name, "_info.TXT", sep = "")

figureInfo <- "3_repdomains.R"

source(paste(sourcePath, "4_figure_info.R", sep = ""))

