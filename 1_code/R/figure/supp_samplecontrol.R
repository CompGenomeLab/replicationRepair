library("ggpubr")
library(reshape2)
library(dplyr)
library("ggplot2")
library("tidyr")
library(grid)
library(gridExtra)
library(patchwork)

# plot layout
layout <- "
AACCDDGG
BBEEFFHH
"

# directories
ts_nts_path <- paste("/home/azgarian/Documents/myprojects/",
                     "replicationRepair/4_output/gitignore/",
                     "2_sample_control/TS_NTS_ratio/", sep = "")
dinuc_path <- paste("~/Documents/myprojects/replicationRepair/4_output/",
                    "gitignore/2_sample_control/dinucleotide_composition/", 
                    sep = "")
corr_path <- paste("~/Documents/myprojects/replicationRepair/4_output/",
                   "gitignore/2_sample_control/correlation/", 
                   sep = "")

# raw data for correlation
rawdata <- read.table(paste(corr_path, "genome_50kb_organized.txt", 
                            sep = ""))

my_data <- dcast(rawdata, V1 + V2 + V3 + V5 + V6 + V7 + V8 + V10 
                 ~ V9, value.var = "V11")

# sample
sample_info <- read.csv(paste("~/Documents/myprojects/replicationRepair/", 
                              "0_data/samples_filtered.csv",  
                              sep = ""))

sample_info_filt <- filter(sample_info, method != "DNA_seq" )

sample_info_agg <- dcast(sample_info_filt, product + release + time 
                         ~ method, value.var = "cell_line")


for(counter in 1:nrow(sample_info_agg)) { 
  
  # correlation plots
  xrcorr_df <- filter(my_data, V5 == sample_info_agg$product[[counter]] & 
                        V6 == "XR_seq" & 
                        V7 == sample_info_agg$release[[counter]] & 
                        V8 == sample_info_agg$time[[counter]])
  
  xrcorrplot <- ggscatter(xrcorr_df, x = "A", y = "B", 
                          add = "reg.line", conf.int = TRUE, 
                          cor.coef = TRUE, cor.method = "spearman",
                          xlab = "Replicate B", 
                          ylab = "XR-seq \n\nReplicate A") +
    theme(axis.title.x = element_text(size = 16, face = "bold"),
          axis.title.y = element_text(size = 16, face = "bold"))
  
  dscorr_df <- filter(my_data, V5 == sample_info_agg$product[[counter]] & 
                        V6 == "Damage_seq" & 
                        V7 == sample_info_agg$release[[counter]] & 
                        V8 == sample_info_agg$time[[counter]])
  
  dscorrplot <- ggscatter(dscorr_df, x = "A", y = "B", 
                          add = "reg.line", conf.int = TRUE, 
                          cor.coef = TRUE, cor.method = "spearman",
                          xlab = "Replicate B", 
                          ylab = "Damage-seq \n\nReplicate A") +
    theme(axis.title.x = element_text(size = 16, face = "bold"),
          axis.title.y = element_text(size = 16, face = "bold"))
  
  # dinucleotide composition
  dinuc_xra <- filter(sample_info_filt, 
                      product == sample_info_agg$product[[counter]] & 
                        method == "XR_seq" & replicate == "A" &
                        release == sample_info_agg$release[[counter]] & 
                        time == sample_info_agg$time[[counter]])
  
  dinuc_xra_df <- read.table(paste(dinuc_path, dinuc_xra$file_name,
                                   "_cutadapt_sorted_26_dinucleotideTable.txt", 
                                   sep = ""), 
                             header = TRUE)
  
  dinuc_xrb <- filter(sample_info_filt, 
                      product == sample_info_agg$product[[counter]] & 
                        method == "XR_seq" & replicate == "B" &
                        release == sample_info_agg$release[[counter]] & 
                        time == sample_info_agg$time[[counter]])
  
  dinuc_xrb_df <- read.table(paste(dinuc_path, dinuc_xrb$file_name,
                                   "_cutadapt_sorted_26_dinucleotideTable.txt", 
                                   sep = ""), 
                             header = TRUE)
  
  dinuc_dsa <- filter(sample_info_filt, 
                      product == sample_info_agg$product[[counter]] & 
                        method == "Damage_seq" & replicate == "A" &
                        release == sample_info_agg$release[[counter]] & 
                        time == sample_info_agg$time[[counter]])
  
  dinuc_dsa_df <- read.table(paste(dinuc_path, dinuc_dsa$file_name,
                                   "_cutadapt_sorted_10_dinucleotideTable.txt", 
                                   sep = ""), 
                             header = TRUE)
  
  dinuc_dsb <- filter(sample_info_filt, 
                      product == sample_info_agg$product[[counter]] & 
                        method == "Damage_seq" & replicate == "B" &
                        release == sample_info_agg$release[[counter]] & 
                        time == sample_info_agg$time[[counter]])
  
  dinuc_dsb_df <- read.table(paste(dinuc_path, dinuc_dsb$file_name,
                                   "_cutadapt_sorted_10_dinucleotideTable.txt", 
                                   sep = ""), 
                             header = TRUE)
  
  dinuc_samp <- list(dinuc_xra_df, dinuc_xrb_df, 
                     dinuc_dsa_df, dinuc_dsb_df)
  
  
  for ( di in 1:length(dinuc_samp) ) {
    
    dinucleotide_table <- dinuc_samp[[di]]
    
    if (di == 1 | di == 3){ myrep <- "A" }
    if (di == 2 | di == 4){ myrep <- "B" }
    
    # rename columns and get their order
    x_order <- c()
    for (i in 2:ncol(dinucleotide_table)) {
      
      colnames(dinucleotide_table)[i] <- c(paste(i - 1, "-", i, sep = ""))
      
      x_order <- c(x_order, paste(i - 1, "-", i, sep = "")) 
      
    }
    
    dt_organized <- dinucleotide_table %>% gather(Position, count, 
                                                  2:ncol(dinucleotide_table))
    
    # reorganize
    colnames(dt_organized) <- c("dinucleotides", "positions", "counts")
    
    dt_organized$freq = 100*dt_organized$counts/sum(dinucleotide_table[2])
    
    dt_organized <- filter(dt_organized, dinucleotides == "CC" |
                             dinucleotides == "CT" | dinucleotides == "TC" |
                             dinucleotides == "TT")
    
    dt_organized$dinucleotides = factor(dt_organized$dinucleotides, 
                                        levels = c("CC", "CT", "TC", "TT"))
    
    dt_organized$positions = factor(dt_organized$positions, 
                                    levels = x_order)
    
    # plot
    
    temp_plotname <- ggplot(dt_organized, aes(x = positions, y = freq, 
                                              fill = dinucleotides)) + 
      geom_bar(stat = "identity") +
      ylim(0,101) +
      scale_fill_manual(values = c("forestgreen", "royalblue3",
                                   "gold1", "mediumvioletred")) +
      ylab("Frequency (%)") + xlab(paste("Replicate", myrep)) + 
      theme_pubr() +
      theme(axis.title.x = element_text(size = 16, face = "bold"),
            axis.title.y = element_text(size = 16, face = "bold"),
            axis.text.x = element_text(size = 14, vjust = 0.6, angle = 65),
            axis.text.y = element_text(size = 14, vjust = 0.1),
            strip.text.x = element_text(size = 14),
            strip.text.y = element_text(size = 14, angle = 360),
            legend.title = element_blank(),
            legend.text = element_text(size = 14),
            plot.caption = element_text(size = 14, 
                                        face = "italic"),
            plot.title = element_text(size = 16, face = "bold"),
            plot.subtitle = element_text(size = 16, face = "bold"))
    
    plotname <- paste("dinuc_plot", di, sep = "")
    
    assign(plotname, temp_plotname)
    
  }
  
  # ts/nts
  ts_nts_xra <- filter(sample_info_filt, 
                       product == sample_info_agg$product[[counter]] & 
                         method == "XR_seq" & replicate == "A" &
                         release == sample_info_agg$release[[counter]] & 
                         time == sample_info_agg$time[[counter]])
  
  ts_nts_xra_df <- read.table(paste(ts_nts_path, ts_nts_xra$file_name,
                                    "_cutadapt_sorted_TSoverNTScount.txt", 
                                    sep = ""), 
                              header = FALSE)
  
  ts_nts_xra_df$replicate <- "A"
  
  ts_nts_xrb <- filter(sample_info_filt, 
                       product == sample_info_agg$product[[counter]] & 
                         method == "XR_seq" & replicate == "B" &
                         release == sample_info_agg$release[[counter]] & 
                         time == sample_info_agg$time[[counter]])
  
  ts_nts_xrb_df <- read.table(paste(ts_nts_path, ts_nts_xrb$file_name,
                                    "_cutadapt_sorted_TSoverNTScount.txt",  
                                    sep = ""), 
                              header = FALSE)
  
  ts_nts_xrb_df$replicate <- "B"
  
  ts_nts_xr_df <- rbind(ts_nts_xra_df, ts_nts_xrb_df)
  
  ts_nts_dsa <- filter(sample_info_filt, 
                       product == sample_info_agg$product[[counter]] & 
                         method == "Damage_seq" & replicate == "A" &
                         release == sample_info_agg$release[[counter]] & 
                         time == sample_info_agg$time[[counter]])
  
  ts_nts_dsa_df <- read.table(paste(ts_nts_path, ts_nts_dsa$file_name,
                                    "_cutadapt_sorted_TSoverNTScount.txt",  
                                    sep = ""), 
                              header = FALSE)
  
  ts_nts_dsa_df$replicate <- "A"
  
  ts_nts_dsb <- filter(sample_info_filt, 
                       product == sample_info_agg$product[[counter]] & 
                         method == "Damage_seq" & replicate == "B" &
                         release == sample_info_agg$release[[counter]] & 
                         time == sample_info_agg$time[[counter]])
  
  ts_nts_dsb_df <- read.table(paste(ts_nts_path, ts_nts_dsb$file_name,
                                    "_cutadapt_sorted_TSoverNTScount.txt",  
                                    sep = ""), 
                              header = FALSE)
  
  ts_nts_dsb_df$replicate <- "B"
  
  ts_nts_ds_df <- rbind(ts_nts_dsa_df, ts_nts_dsb_df)
  
  ts_nts_samp <- list(ts_nts_xr_df, ts_nts_ds_df)
  
  for ( ts in 1:length(ts_nts_samp) ) {
    
    ts_nts_df <- ts_nts_samp[[ts]]
    
    colnames(ts_nts_df) <- c("chromosomes", "s_point", "e_point", "gene_id", 
                             "strand", "TS", "NTS","replicate") 
    
    ts_nts_df$TSoverNTS <- ts_nts_df$TS / ts_nts_df$NTS
    
    #### plot ####
    
    temp_ts_plotname <- ggplot(ts_nts_df, aes(x = replicate, 
                                              y = log2(TSoverNTS))) + 
      geom_boxplot(outlier.shape = NA) + 
      xlab("Replicate") + ylab("log2 normalized TS/NTS") +
      ylim(-1, 1) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "red") 
    
    #### theme ####
    
    temp_ts_plotname <- temp_ts_plotname + theme_pubr() +
      theme(axis.title.x = element_text(size = 16, face = "bold"),
            axis.title.y = element_text(size = 16, face = "bold"),
            axis.text.x = element_text(size = 16),
            axis.text.y = element_text(size = 14, vjust = 0.1),
            strip.text.x = element_text(size = 14),
            strip.text.y = element_text(size = 14, angle = 360),
            legend.title = element_text(size = 14, face = "bold"),
            legend.text = element_text(size = 14),
            plot.caption = element_text(size = 14, 
                                        face = "italic"),
            plot.title = element_text(size = 16, face = "bold"),
            plot.subtitle = element_text(size = 16, face = "bold"))
    
    plotname <- paste("ts_nts_plot", ts, sep = "")
    
    assign(plotname, temp_ts_plotname)
    
  }
  
  
  # combine and save plots
  
  xrcorrplot + dscorrplot + dinuc_plot1 + dinuc_plot2 + dinuc_plot3 + 
    dinuc_plot4 + ts_nts_plot1 + ts_nts_plot2 +
    plot_layout(design = layout, guides = "collect", tag_level = 'new') & 
    theme(legend.position = 'bottom') 
  
  ggsave(paste("~/Desktop/sample_control_", sample_info_agg$product[[counter]], 
               "_", sample_info_agg$release[[counter]], "_", 
               sample_info_agg$time[[counter]],
               ".png", sep = ""), width = 597, height = 310, units = "mm")
  
}

