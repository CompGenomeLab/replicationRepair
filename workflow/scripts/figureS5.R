#### Packages and Libraries ####

library(stringr)
library(ggpubr)
library(ggplot2)
library(dplyr)
library(reshape2)
library(patchwork)
library(grid)
library(argparser)
library(scales)
library(zoo)
set.seed(1) 

######## Arguments ##########
p <- arg_parser("producing the supplementary figure 5")
p <- add_argument(p, "--df", help="region file with read counts")
p <- add_argument(p, "--df_sim", help="region file with simulated read counts")
p <- add_argument(p, "--data_prefix", help="name prefix of the dataframes that generate the plots")
p <- add_argument(p, "-o", help="output")

# Parse the command line arguments
argv <- parse_args(p)

sample <- argv$df

sample_sim <- argv$df_sim

#### Default Plot Format ####

source("workflow/scripts/plot_format.R")

options(scipen=999)


#### Functions ####

source("workflow/scripts/functions.R")


#### Main ####

# for plot A
sample_df <- read.delim( sample, header = F )
colnames(sample_df) <- c("chromosomes", "start_position", "end_position", 
                            "dataset", "score", "dataset_strand", "counts", 
                            "sample_names", "file_names", "layout", "cell_line", 
                            "product", "method", "uv_exposure", "treatment", 
                            "phase", "time_after_exposure", "replicate", 
                            "project", "sample_source", "sample_strand", 
                            "mapped_reads", "RPKM")

df_rr <- repair_rate( sample_df )

sample_df_sim <- read.delim( sample_sim, header = F )

colnames(sample_df_sim) <- c("chromosomes", "start_position", "end_position", 
                             "dataset", "score", "dataset_strand", "counts", 
                             "sample_names", "file_names", "layout", "cell_line", 
                             "product", "method", "uv_exposure", "treatment", 
                             "phase", "time_after_exposure", "replicate", 
                             "project", "sample_source", "sample_strand", 
                             "mapped_reads", "RPKM")

df_rr_sim <- repair_rate( sample_df_sim )

colnames(df_rr)[12] <- "real"
colnames(df_rr_sim)[12] <- "sim"

df_rr_rs <- merge(df_rr, df_rr_sim, by = c("chromosomes", "start_position", 
                                           "end_position", "dataset", "score",
                                           "dataset_strand", "product", 
                                           "phase", "time_after_exposure",
                                           "sample_strand", "replicate" ))
df_rr_rs$real_sim <- df_rr_rs$real / df_rr_rs$sim

df_rr_rs_score <- data.frame()
for (chr in unique(df_rr_rs$chromosomes)){
  chr_data <- filter(df_rr_rs, chromosomes == chr)
  chr_data <- chr_data[order(chr_data$start_position),]
  
  zoo_mean<-aggregate(chr_data$score, list(chr_data$start_position), mean)
  #Make zoo object of data
  temp.zoo<-zoo(zoo_mean$x,zoo_mean$Group.1)
  
  #Calculate moving average with window 3 and make first and last value as NA (to ensure identical length of vectors)
  score_sma<-rollmean(temp.zoo, 10, na.rm = TRUE)
  
  score_sma <- data.frame(score_sma)
  score_sma$start_position <- rownames(score_sma)
  score_sma$chromosomes <- chr
  
  chr_data_merged <- merge(chr_data, score_sma, by=c("chromosomes", "start_position"))
  
  df_rr_rs_score <- rbind(df_rr_rs_score, chr_data_merged)
}

df_rr_rs_dcast_plus_minus <- 
  dcast(df_rr_rs_score, dataset + chromosomes + start_position + end_position + 
          score + score_sma + dataset_strand + product + phase + 
          time_after_exposure + replicate ~ sample_strand, 
        value.var = "real_sim")
df_rr_rs_dcast_plus_minus$plus_minus <- df_rr_rs_dcast_plus_minus$`+` / df_rr_rs_dcast_plus_minus$`-`

df_rr_rs_dcast_ear_la <- 
  dcast(df_rr_rs_score, dataset + chromosomes + start_position + end_position + 
          score + score_sma + dataset_strand + product + sample_strand + 
          time_after_exposure + replicate ~ phase, value.var = "real_sim")
df_rr_rs_dcast_ear_la$ear_la <- df_rr_rs_dcast_ear_la$early / df_rr_rs_dcast_ear_la$late

# filtering samples 
pA_data <- filter(df_rr_rs_score, 
              replicate == "_", time_after_exposure == "12",
              phase != "async", product == "CPD", chromosomes == "chr1" | 
              chromosomes == "chr19", start_position < 50000000 )

pB_data <- filter(df_rr_rs_dcast_ear_la, 
                replicate == "_", time_after_exposure == "12",
                product == "CPD", chromosomes == "chr1" | 
                chromosomes == "chr19", start_position < 50000000)


#### Plot P.A ####

write.table(pA_data, file = paste0(argv$data_prefix, "A.csv"), 
            quote = FALSE, row.names = FALSE, sep = ",")

# create the plot
p.A <- ggplot(pA_data, aes(x = start_position, y = score_sma, color = log2(real_sim))) + 
  geom_line(size=3) +
  #geom_ma(ma_fun = SMA, n = 30, data=data, na.rm = TRUE) +
  #geom_smooth(se=FALSE) +
  facet_grid(~phase~chromosomes,
             labeller = labeller(product = product_labs, 
                                 time_after_exposure = taex_labs,
                                 phase = phase_labs)) +
  scale_color_gradientn( colours = c("red3","white","blue4"),
                       values=rescale(c(-4, 0,2)), limits=c(-4,2)) +
  ylab("Replication Timing (Early/Late)") +
  xlab("Chromosome Position (Mb)") +
  scale_x_continuous(breaks = c(0, 25000000, 50000000),
                     labels = c("0","25", "50")) +
  labs(color="n. Repair\nRate (log2)")

# adding and overriding the default plot format
p.A <- p.A + p_format + theme(legend.position = "right",
                              axis.title.x = element_blank())


#### Plot P.B ####

write.table(pB_data, file = paste0(argv$data_prefix, "B.csv"), 
            quote = FALSE, row.names = FALSE, sep = ",")

# create the plot
p.B <- ggplot(pB_data, aes(x = start_position, y = score_sma, color = log2(ear_la))) + 
  geom_line(size=3) +
  #geom_ma(ma_fun = SMA, n = 30, data=data, na.rm = TRUE) +
  #geom_smooth(se=FALSE) +
  facet_grid(~chromosomes,
             labeller = labeller(product = product_labs, 
                                 time_after_exposure = taex_labs,
                                 phase = phase_labs)) +
  scale_color_gradientn( colours = c("red3","white","blue4"),
                         values=rescale(c(-1, 0,1)), limits=c(-1,1)) +
  ylab("Replication Timing (Early/Late)") +
  xlab("Chromosome Position (Mb)") +
  scale_x_continuous(breaks = c(0, 25000000, 50000000),
                     labels = c("0","25", "50")) +
  labs(color="Early/Late \nS Phase \nn. Repair\nRate (log2)")

# adding and overriding the default plot format
p.B <- p.B + p_format + theme(legend.position = "right") + 
  theme(strip.text.x = element_blank())

layout <- 
"
A
A
B
"

p.A / p.B + plot_layout(design = layout) &
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size = 12, face="bold"))

ggsave(argv$o, width = 22, height = 18, units = "cm") 
