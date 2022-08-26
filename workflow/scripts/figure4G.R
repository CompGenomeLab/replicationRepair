#### Packages and Libraries ####

library(stringr)
library(ggpubr)
library(ggplot2)
library(dplyr)
library(reshape2)
library(patchwork)
library(grid)
library(argparser)

######## Arguments ##########
p <- arg_parser("producing the figure 1")
p <- add_argument(p, "--df", help="regions with XR and DS values")
p <- add_argument(p, "--counts", help="nucleotide content of regions")
p <- add_argument(p, "--data_prefix", help="name prefix of the dataframes that generate the plots")
p <- add_argument(p, "-o", help="output")
p <- add_argument(p, "--log", help="log file")

# Parse the command line arguments
argv <- parse_args(p)

sample <- argv$df
counts <- argv$counts

#### Default Plot Format ####

source("workflow/scripts/plot_format.R")

sample_strand_labs <- c("Reverse Strand", "Forward Strand")
names(sample_strand_labs) <- c("-", "+")

quant_labs <- c("A-rich Regions", "T-rich Regions")
names(quant_labs) <- c(1, 4)

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

df_rr <- dcast(sample_df, chromosomes + start_position + end_position + 
                 product + sample_strand + time_after_exposure + phase +
                 replicate ~ method, value.var = "RPKM")

df_rr$xr_ds <- df_rr$XR_seq / df_rr$Damage_seq

sample_counts <- read.delim( counts, header = T )

#sample_counts = sample_counts[,c(1,2,3,4,5,6,7)]

df_rr_counts <- merge(df_rr, sample_counts, by.x = c("chromosomes", 
                                                        "start_position", 
                                                        "end_position"),
                      by.y = c ("chr","start","end"))

df_rr_counts$relativeTs <- log2(df_rr_counts$t / df_rr_counts$a)
df_rr_counts$relativeCs <- log2(df_rr_counts$c / df_rr_counts$g)

df_rr_counts <- filter(df_rr_counts, relativeTs != "Na", relativeTs != "NaN",
                       relativeCs != "Na", relativeCs != "NaN")

df_rr_counts <- within(df_rr_counts, relT_quant <- 
                         as.integer(cut(relativeTs, quantile(relativeTs, 
                                                             probs=0:4/4), 
                                        include.lowest=TRUE)))

df_rr_counts <- within(df_rr_counts, relC_quant <- 
                         as.integer(cut(relativeCs, quantile(relativeCs, 
                                                             probs=0:4/4), 
                                        include.lowest=TRUE)))

df_rr_counts$t_quant <- 0
df_rr_counts$t_quant[df_rr_counts$t > 35] <- 4
df_rr_counts$t_quant[30 < df_rr_counts$t & df_rr_counts$t <= 35] <- 3
df_rr_counts$t_quant[25 < df_rr_counts$t &  df_rr_counts$t <= 30] <- 2
df_rr_counts$t_quant[df_rr_counts$t <= 25] <- 1

df_rr_counts$c_quant <- 0
df_rr_counts$c_quant[df_rr_counts$c > 35] <- 4
df_rr_counts$c_quant[30 < df_rr_counts$c & df_rr_counts$c <= 35] <- 3
df_rr_counts$c_quant[25 < df_rr_counts$c &  df_rr_counts$c <= 30] <- 2
df_rr_counts$c_quant[df_rr_counts$c <= 25] <- 1

df_rr_counts$g_quant <- 0
df_rr_counts$g_quant[df_rr_counts$g > 35] <- 4
df_rr_counts$g_quant[30 < df_rr_counts$g & df_rr_counts$g <= 35] <- 3
df_rr_counts$g_quant[25 < df_rr_counts$g &  df_rr_counts$g <= 30] <- 2
df_rr_counts$g_quant[df_rr_counts$g <= 25] <- 1

df_rr_counts$a_quant <- 0
df_rr_counts$a_quant[df_rr_counts$a > 35] <- 4
df_rr_counts$a_quant[30 < df_rr_counts$a & df_rr_counts$a <= 35] <- 3
df_rr_counts$a_quant[25 < df_rr_counts$a &  df_rr_counts$a <= 30] <- 2
df_rr_counts$a_quant[df_rr_counts$a <= 25] <- 1

# filtering samples 
data <- filter(df_rr_counts, replicate == "_", n < 5, 
               time_after_exposure == "12", product == "CPD", phase == "async")

data$sample_strand <- factor(
  data$sample_strand, levels = c("+","-"))

data_relT <- filter(data, relT_quant == 1 |  relT_quant == 4)

data_relC <- filter(data, relC_quant == 1 |  relC_quant == 4)

#### relative T #####

write.table(data_relT, file = paste0(argv$data_prefix, ".csv"), quote = FALSE, 
            row.names = FALSE, sep = ",")

p.A <- ggplot(data_relT, aes(x = sample_strand, y = Damage_seq, 
                        fill=sample_strand)) +
  geom_boxplot(outlier.shape = NA) + 
  facet_grid(~relT_quant, labeller = 
               labeller(relT_quant = quant_labs)) +
  ylab("Damage-seq RPKM") +
  scale_fill_manual(values=c("white", "white")) +
  scale_x_discrete(labels = c("Forward\nStrand", "Reverse\nStrand")) + 
  guides(fill = "none") +
  ylim(0,0.75)

p.A <- p.A + p_format + 
  theme(strip.background = element_blank(),
        axis.title.x=element_blank()) +
  stat_compare_means(comparisons = list( c(1,2)),
                     label.y = .4, label = "p.signif")

p.B <- ggplot(data_relT, aes(x = sample_strand, y = XR_seq, 
                        fill=sample_strand)) +
  geom_boxplot(outlier.shape = NA) + 
  facet_grid(~relT_quant, labeller = 
               labeller(relT_quant = quant_labs)) +
  ylab("XR-seq RPKM") +
  scale_fill_manual(values=c("white", "white")) +
  scale_x_discrete(labels = c("Forward\nStrand", "Reverse\nStrand")) + 
  guides(fill = "none") +
  ylim(0,0.75)

p.B <- p.B + p_format + 
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.title.x=element_blank()) +
  stat_compare_means(comparisons = list( c(1,2)),
                     label.y = .5, label = "p.signif")

p.C <- ggplot(data_relT, aes(x = sample_strand, y = log2(xr_ds), 
                        fill=sample_strand)) +
  geom_boxplot(outlier.shape = NA) + 
  facet_grid(~relT_quant, labeller = 
               labeller(relT_quant = quant_labs)) +
  xlab("") + 
  ylab("Repair Rate (log2)") +
  scale_fill_manual(values=c("white", "white")) +
  scale_x_discrete(labels = c("Forward\nStrand", "Reverse\nStrand")) + 
  guides(fill = "none") +
  ylim(-2.5,4)

p.C <- p.C + p_format + 
  theme(strip.text.x = element_blank()) +
  stat_compare_means(comparisons = list( c(1,2)),
                     label.y = 3, label = "p.signif")

layout <- "
AAAA
BBBB
CCCC
"

p.A + p.B + p.C +
  plot_layout(design = layout) & 
  theme(plot.tag = element_text(hjust = -0.2, size = 10, face="bold"))

ggsave( argv$o, width = 10, height = 15, units = "cm" )
