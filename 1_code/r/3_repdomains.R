#### set paths and variables ####

# name of the data file and its run date  
#fr_name <- "repdomains"
#date <- "[2020.02.10]"

library(ggpubr) # for stat_compare_means()

# path of the codes
sourcePath <- "~/Documents/myprojects/replicationRepair/1_code/R/" 

# directory of figures
figurePath <- paste("~/Documents/My_Projects/Project_Repair_Replication/", 
                    "Results/Repdomains/", sep = "")

options(scipen = 999) # turn off scientific notation

# file name of data information
dataInfo <- paste(date, "final_report_", fr_name, "_info.TXT", sep = "") 

# name of the script
figureInfo <- "3_repdomains.R" 

# pre-analysis of data file
source(paste(sourcePath, "2_report_sub_dfs.R", sep = ""))

# functions
source(paste(sourcePath, "4_functions.R", sep = ""))

# add plot format
source(paste(sourcePath, "4_plot_format.R", sep = ""))

# date of today
dateout <- output_date()


#### Repair Rate of Chromosomes ####

# rearrange
df <- fr_xr_ds
df$chromosomes <- factor(df$chromosomes, levels = c("chr1", "chr2", "chr3", 
                                                  "chr4", "chr5", "chr6", 
                                                  "chr7", "chr8", "chr9", 
                                                  "chr10", "chr11", "chr12", 
                                                  "chr13", "chr14", "chr15", 
                                                  "chr16", "chr17", "chr18", 
                                                  "chr19", "chr20", "chr21", 
                                                  "chr22", "chrX"))

# filter the data
d <- filter(df, dataset == "ERD" | dataset == "LRD", replicate == "A", 
            product == "CPD", time_after_exposure == "12")

# create the plot 
p <- ggplot(d, aes(x = chromosomes, y = log2(xr_ds))) + 
  geom_boxplot(aes(fill = dataset), outlier.shape = NA) +
  geom_hline(yintercept = -1, linetype="dashed", color = "darkred") +
  geom_hline(yintercept = 1, linetype="dashed", color = "darkgreen") +
  facet_grid(phase~time_after_exposure~replicate~product,
             labeller = labeller(product = product_labs, 
                                 time_after_exposure = taex_labs, 
                                 replicate = rep_labs, phase = phase_labs)) +
  xlab("Chromosomes") + ylab(fr_xr_ds_lab) +
  scale_fill_manual(name="", values = repdomain_colors) + ylim(-3,2)
p <- p +  p_format + theme(axis.text.x = 
                             element_text(size = 12, vjust = 0.6, angle = 60))
# p # visualize

# save figure
figureName <- paste(dateout, "Repair_over_Damage_Value_of_", 
                    fr_name, "_for_Every_Chromosome_CPD_12.pdf", sep = "")

fig_save( figurePath, figureName )

# save figure info
source(paste(sourcePath, "4_figure_info.R", sep = ""))


#### Repair Rate of ERD and LRD ####

# rearrange
df <- fr_xr_ds
phase_labs <- c("Async.", "Early Phase", "Late Phase")
names(phase_labs) <- c("async", "early", "late")

#### filter the data
d = filter(df, dataset == "ERD" | dataset == "LRD", time_after_exposure == "12")

# create the plot
p <- ggplot(d, aes(x = dataset, y = log2(xr_ds))) + 
  geom_boxplot(aes(fill = dataset), outlier.shape = NA) +
  facet_grid(product~time_after_exposure~replicate~phase, 
             labeller = labeller(product = product_labs, 
                                 time_after_exposure = taex_labs, 
                                 phase = phase_labs, replicate = rep_labs)) + 
  xlab(repdo_lab) + ylab(fr_xr_ds_lab) +
  geom_hline(yintercept = 0, linetype="dashed", color="red") +
  scale_fill_manual(name="", values = repdomain_colors) +
  ylim(-2,2) 

p <- p + p_format + stat_compare_means(method = "t.test", 
                                       label.x = 1.3, label.y = 1.8)
p # visualize

# save figure
figureName <- paste(dateout, "Repair_over_Damage_Value_of_", 
                    fr_name, "_for_Every_Sample.pdf", sep = "")

fig_save( figurePath, figureName )

# save figure info
source(paste(sourcePath, "4_figure_info.R", sep = ""))
