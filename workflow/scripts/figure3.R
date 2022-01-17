#### Packages and Libraries ####

library(ggpubr)
library(ggplot2)
library(stringr)
library(reshape2)
library(dplyr)
library(patchwork)
library(argparser)

######## Arguments ##########
p <- arg_parser("producing the figure 3")
p <- add_argument(p, "--df", help="region file with read counts")
p <- add_argument(p, "--df_sim", help="region file with simulated read counts")
p <- add_argument(p, "--dtype", help="damage type of the sample to visualize")
p <- add_argument(p, "--o1", help="output")
p <- add_argument(p, "--o2", help="output (paired t-test)")

# Parse the command line arguments
argv <- parse_args(p)

sample_csv <- argv$df

sample_csv_sim <- argv$df_sim

#### Variables ####

rep <- "_" # replicate
dprod <- argv$dtype # damage product (CPD/64_PP)

#sample_csv <- paste("/Users/azgarian/Desktop/replication_final/",
#                    "final_reports_hg19_chromhmm_repdomains_hela_windows_chr.txt", 
#                    sep = "")

#sample_csv_sim <- "/Users/azgarian/Desktop/replication_final/final_reports_sim_hg19_chromhmm_repdomains_hela_windows_chr.txt"

#### Default Plot Format ####

source("workflow/scripts/plot_format.R")


#### Fuctions ####

source("workflow/scripts/functions.R")

#### Main ####

# for both plot A and plot B
sample_df <- read.delim( sample_csv, header = F )

colnames(sample_df) <- c("chromosomes", "start_position", "end_position", 
                            "dataset", "score", "dataset_strand", "counts", 
                            "sample_names", "file_names", "layout", "cell_line", 
                            "product", "method", "uv_exposure", "treatment", 
                            "phase", "time_after_exposure", "replicate", 
                            "project", "sample_source", "sample_strand", 
                            "mapped_reads", "RPKM")

sample_df_filt <- filter(sample_df, product == dprod, replicate == rep, 
                         phase != "async", time_after_exposure == "12")

sample_df_rmlow <- rmv_low_counted_regions(sample_df_filt, 2)

df_rr <- repair_rate( sample_df_rmlow )

colnames(df_rr)[12] <- "real"

sample_df_sim <- read.delim( sample_csv_sim, header = F )

colnames(sample_df_sim) <- c("chromosomes", "start_position", "end_position", 
                         "dataset", "score", "dataset_strand", "counts", 
                         "sample_names", "file_names", "layout", "cell_line", 
                         "product", "method", "uv_exposure", "treatment", 
                         "phase", "time_after_exposure", "replicate", 
                         "project", "sample_source", "sample_strand", 
                         "mapped_reads", "RPKM")

sample_df_sim_filt <- filter(sample_df_sim, product == dprod, replicate == rep, 
                         phase != "async", time_after_exposure == "12")

sample_df_sim_rmlow <- rmv_low_counted_regions(sample_df_sim_filt, 2)

df_rr_sim <- repair_rate( sample_df_sim_rmlow )

colnames(df_rr_sim)[12] <- "sim"

df_rr_rs <- merge(df_rr, df_rr_sim, by = c("chromosomes", "start_position", 
                                            "end_position", "dataset", "score",
                                            "dataset_strand", "product", 
                                            "phase", "time_after_exposure", 
                                            "sample_strand", "replicate" ))
df_rr_rs$xr_ds <- df_rr_rs$real / df_rr_rs$sim

# for plot A
pA_df_rr <- chrState_naming( df_rr_rs, chr_states, general_states, 
                             chrState2generalState )

pA_df_rr$state_short <- NA
for (i in 1:length(chr_states)) {
  pA_df_rr <- within(pA_df_rr, state_short[states == general_states[i]] <- 
                       general_states_short[i])
}

# for plot B
df_rr_ear_la <- rrEarly_rrLate( df_rr_rs )
pB_df_rr_ear_la <- chrState_naming( df_rr_ear_la, chr_states, general_states, 
                                    chrState2generalState )

pB_df_rr_ear_la$state_short <- NA
for (i in 1:length(chr_states)) {
  pB_df_rr_ear_la <- within(pB_df_rr_ear_la, state_short[states == general_states[i]] <- 
                       general_states_short[i])
}

#### Filtering Samples ####

# for plot A
pA_filt <- filter(pA_df_rr, xr_ds != "NaN", xr_ds != "Inf", 
                  dataset_strand != "DTZ", dataset_strand != "UTZ", 
                  product == dprod, replicate == rep, phase != "async", 
                  time_after_exposure == "12")

pA_filt$state_short <- factor(pA_filt$state_short, levels=general_states_short)

pA_data <- aggregate(x = pA_filt$xr_ds, by = list(pA_filt$dataset, 
                                                  pA_filt$dataset_strand, 
                                                  pA_filt$phase, 
                                                  pA_filt$states, 
                                                  pA_filt$state_short,
                                                  pA_filt$chromosomes), 
                     FUN = "mean")

pA_data$Group.5 <- factor(pA_data$Group.5, levels=general_states_short)

# for plot B
pB_filt <- filter(pB_df_rr_ear_la, ear_la != "NaN", ear_la != "Inf", 
                  dataset_strand != "DTZ", dataset_strand != "UTZ", 
                  product == dprod, replicate == rep,
                  time_after_exposure == "12")
pB_filt$state_short <- factor(pB_filt$state_short, levels=general_states_short)

pB_data <- aggregate(x = pB_filt$ear_la, by = list(pB_filt$dataset, 
                                                   pB_filt$dataset_strand, 
                                                   pB_filt$states, 
                                                   pB_filt$state_short,
                                                   pB_filt$chromosomes), 
                     FUN = "mean")

pB_data$Group.4 <- factor(pB_data$Group.4, levels=general_states_short)

#### Plot A ####

pA_filt$phase[pA_filt$phase == "early"] <- "Early S Phase"
pA_filt$phase[pA_filt$phase == "late"] <- "Late S Phase"

# create the plot
p.A <- ggplot(pA_filt, aes(x = state_short, y = log2(xr_ds))) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_boxplot(aes(fill=factor(states), linetype = factor(phase)), 
               outlier.shape = NA, lwd=0.5) +
  facet_grid(dataset_strand~.) +
  xlab("") + ylab("n. Repair Rate (RR) (log2)") +
  scale_fill_manual(name = "Chromatin States", values = state_colors) +
  scale_linetype_manual(name = "Phases", values = c("Early S Phase" =  "solid", 
                                                    "Late S Phase" = "dashed")) +
  ylim(-4,5) +
  guides(size = "none") 

# adding and overriding the default plot format
p.A <- p.A + p_format +
  theme(panel.border=element_rect(size=1, fill = NA),
        strip.text = element_blank(),
        legend.position = "right", 
        axis.text.x = element_text(angle = 60, vjust = 0.95, hjust = 0.95))

ggsave(argv$o1, width = 22, height = 18, units = "cm")

#### Plot A_ttest ####

# create the plot
p.A_ttest <- ggplot(pA_filt, aes(x = phase, y = log2(xr_ds))) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_boxplot(aes(fill=factor(states), linetype = factor(phase)), 
               outlier.shape = NA, lwd=0.5) +
  facet_grid(~dataset_strand~state_short) +
  xlab("") + ylab("n. Repair Rate (RR) (log2)") +
  scale_fill_manual(name = "Chromatin States", values = state_colors) +
  #ylim(-2,3) +
  guides(size = "none", fill="none") 

# adding and overriding the default plot format
p.A_ttest <- p.A_ttest + p_format + 
  theme(panel.border=element_rect(size=1, fill = NA),
        legend.position = "right", 
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())

p.A_ttest <- p.A_ttest +
  stat_compare_means(method = "t.test", label = "p.signif", paired = T, 
                     comparisons = list(c("Early S Phase", "Late S Phase")), label.y = 4)

ggsave(argv$o2, width = 22, height = 18, units = "cm")
